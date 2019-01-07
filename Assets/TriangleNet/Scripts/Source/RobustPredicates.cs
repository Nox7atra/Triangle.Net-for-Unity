// -----------------------------------------------------------------------
// <copyright file="RobustPredicates.cs">
// Original Triangle code by Jonathan Richard Shewchuk, http://www.cs.cmu.edu/~quake/triangle.html
// Triangle.NET code by Christian Woltering, http://triangle.codeplex.com/
// </copyright>
// -----------------------------------------------------------------------

namespace TriangleNet
{
    using System;
    using TriangleNet.Geometry;
    using TriangleNet.Tools;

    /// <summary>
    /// Adaptive exact arithmetic geometric predicates.
    /// </summary>
    /// <remarks>
    /// The adaptive exact arithmetic geometric predicates implemented herein are described in
    /// detail in the paper "Adaptive Precision Floating-Point Arithmetic and Fast Robust
    /// Geometric Predicates." by Jonathan Richard Shewchuk, see
    /// http://www.cs.cmu.edu/~quake/robust.html
    /// 
    /// The macros of the original C code were automatically expanded using the Visual Studio
    /// command prompt with the command "CL /P /C EXACT.C", see
    /// http://msdn.microsoft.com/en-us/library/8z9z0bx6.aspx
    /// </remarks>
    public class RobustPredicates : IPredicates
    {
        #region Default predicates instance (Singleton)

        private static readonly object creationLock = new object();
        private static RobustPredicates _default;

        /// <summary>
        /// Gets the default configuration instance.
        /// </summary>
        public static RobustPredicates Default
        {
            get
            {
                if (_default == null)
                {
                    lock (creationLock)
                    {
                        if (_default == null)
                        {
                            _default = new RobustPredicates();
                        }
                    }
                }

                return _default;
            }
        }

        #endregion

        #region Static initialization

        private static float epsilon, splitter, resulterrbound;
        private static float ccwerrboundA, ccwerrboundB, ccwerrboundC;
        private static float iccerrboundA, iccerrboundB, iccerrboundC;
        //private static float o3derrboundA, o3derrboundB, o3derrboundC;

        /// <summary>
        /// Initialize the variables used for exact arithmetic.  
        /// </summary>
        /// <remarks>
        /// 'epsilon' is the largest power of two such that 1.0 + epsilon = 1.0 in
        /// floating-point arithmetic. 'epsilon' bounds the relative roundoff
        /// error. It is used for floating-point error analysis.
        ///
        /// 'splitter' is used to split floating-point numbers into two half-
        /// length significands for exact multiplication.
        ///
        /// I imagine that a highly optimizing compiler might be too smart for its
        /// own good, and somehow cause this routine to fail, if it pretends that
        /// floating-point arithmetic is too much like float arithmetic.
        ///
        /// Don't change this routine unless you fully understand it.
        /// </remarks>
        static RobustPredicates()
        {
            float half;
            float check, lastcheck;
            bool every_other;

            every_other = true;
            half = 0.5f;
            epsilon = 1.0f;
            splitter = 1.0f;
            check = 1.0f;
            // Repeatedly divide 'epsilon' by two until it is too small to add to
            // one without causing roundoff.  (Also check if the sum is equal to
            // the previous sum, for machines that round up instead of using exact
            // rounding.  Not that these routines will work on such machines.)
            do
            {
                lastcheck = check;
                epsilon *= half;
                if (every_other)
                {
                    splitter *= 2.0f;
                }
                every_other = !every_other;
                check = 1.0f + epsilon;
            } while ((check != 1.0f) && (check != lastcheck));
            splitter += 1.0f;
            // Error bounds for orientation and incircle tests. 
            resulterrbound = (3.0f + 8.0f * epsilon) * epsilon;
            ccwerrboundA = (3.0f + 16.0f * epsilon) * epsilon;
            ccwerrboundB = (2.0f + 12.0f * epsilon) * epsilon;
            ccwerrboundC = (9.0f + 64.0f * epsilon) * epsilon * epsilon;
            iccerrboundA = (10.0f + 96.0f * epsilon) * epsilon;
            iccerrboundB = (4.0f + 48.0f * epsilon) * epsilon;
            iccerrboundC = (44.0f + 576.0f * epsilon) * epsilon * epsilon;
            //o3derrboundA = (7.0 + 56.0 * epsilon) * epsilon;
            //o3derrboundB = (3.0 + 28.0 * epsilon) * epsilon;
            //o3derrboundC = (26.0 + 288.0 * epsilon) * epsilon * epsilon;
        }

        #endregion

        public RobustPredicates()
        {
            AllocateWorkspace();
        }

        /// <summary>
        /// Check, if the three points appear in counterclockwise order. The result is 
        /// also a rough approximation of twice the signed area of the triangle defined 
        /// by the three points.
        /// </summary>
        /// <param name="pa">Point a.</param>
        /// <param name="pb">Point b.</param>
        /// <param name="pc">Point c.</param>
        /// <returns>Return a positive value if the points pa, pb, and pc occur in 
        /// counterclockwise order; a negative value if they occur in clockwise order; 
        /// and zero if they are collinear.</returns>
        public float CounterClockwise(Point pa, Point pb, Point pc)
        {
            float detleft, detright, det;
            float detsum, errbound;

            Statistic.CounterClockwiseCount++;

            detleft = (pa.x - pc.x) * (pb.y - pc.y);
            detright = (pa.y - pc.y) * (pb.x - pc.x);
            det = detleft - detright;

            if (Behavior.NoExact)
            {
                return det;
            }

            if (detleft > 0.0f)
            {
                if (detright <= 0.0f)
                {
                    return det;
                }
                else
                {
                    detsum = detleft + detright;
                }
            }
            else if (detleft < 0.0f)
            {
                if (detright >= 0.0f)
                {
                    return det;
                }
                else
                {
                    detsum = -detleft - detright;
                }
            }
            else
            {
                return det;
            }

            errbound = ccwerrboundA * detsum;
            if ((det >= errbound) || (-det >= errbound))
            {
                return det;
            }

            Statistic.CounterClockwiseAdaptCount++;
            return CounterClockwiseAdapt(pa, pb, pc, detsum);
        }

        /// <summary>
        /// Check if the point pd lies inside the circle passing through pa, pb, and pc. The 
        /// points pa, pb, and pc must be in counterclockwise order, or the sign of the result 
        /// will be reversed.
        /// </summary>
        /// <param name="pa">Point a.</param>
        /// <param name="pb">Point b.</param>
        /// <param name="pc">Point c.</param>
        /// <param name="pd">Point d.</param>
        /// <returns>Return a positive value if the point pd lies inside the circle passing through 
        /// pa, pb, and pc; a negative value if it lies outside; and zero if the four points 
        /// are cocircular.</returns>
        public float InCircle(Point pa, Point pb, Point pc, Point pd)
        {
            float adx, bdx, cdx, ady, bdy, cdy;
            float bdxcdy, cdxbdy, cdxady, adxcdy, adxbdy, bdxady;
            float alift, blift, clift;
            float det;
            float permanent, errbound;

            Statistic.InCircleCount++;

            adx = pa.x - pd.x;
            bdx = pb.x - pd.x;
            cdx = pc.x - pd.x;
            ady = pa.y - pd.y;
            bdy = pb.y - pd.y;
            cdy = pc.y - pd.y;

            bdxcdy = bdx * cdy;
            cdxbdy = cdx * bdy;
            alift = adx * adx + ady * ady;

            cdxady = cdx * ady;
            adxcdy = adx * cdy;
            blift = bdx * bdx + bdy * bdy;

            adxbdy = adx * bdy;
            bdxady = bdx * ady;
            clift = cdx * cdx + cdy * cdy;

            det = alift * (bdxcdy - cdxbdy)
                + blift * (cdxady - adxcdy)
                + clift * (adxbdy - bdxady);

            if (Behavior.NoExact)
            {
                return det;
            }

            permanent = (Math.Abs(bdxcdy) + Math.Abs(cdxbdy)) * alift
                      + (Math.Abs(cdxady) + Math.Abs(adxcdy)) * blift
                      + (Math.Abs(adxbdy) + Math.Abs(bdxady)) * clift;
            errbound = iccerrboundA * permanent;
            if ((det > errbound) || (-det > errbound))
            {
                return det;
            }

            Statistic.InCircleAdaptCount++;
            return InCircleAdapt(pa, pb, pc, pd, permanent);
        }

        /// <summary>
        /// Return a positive value if the point pd is incompatible with the circle 
        /// or plane passing through pa, pb, and pc (meaning that pd is inside the 
        /// circle or below the plane); a negative value if it is compatible; and 
        /// zero if the four points are cocircular/coplanar. The points pa, pb, and 
        /// pc must be in counterclockwise order, or the sign of the result will be 
        /// reversed.
        /// </summary>
        /// <param name="pa">Point a.</param>
        /// <param name="pb">Point b.</param>
        /// <param name="pc">Point c.</param>
        /// <param name="pd">Point d.</param>
        /// <returns>Return a positive value if the point pd lies inside the circle passing through 
        /// pa, pb, and pc; a negative value if it lies outside; and zero if the four points 
        /// are cocircular.</returns>
        public float NonRegular(Point pa, Point pb, Point pc, Point pd)
        {
            return InCircle(pa, pb, pc, pd);
        }

        /// <summary>
        /// Find the circumcenter of a triangle.
        /// </summary>
        /// <param name="org">Triangle point.</param>
        /// <param name="dest">Triangle point.</param>
        /// <param name="apex">Triangle point.</param>
        /// <param name="xi">Relative coordinate of new location.</param>
        /// <param name="eta">Relative coordinate of new location.</param>
        /// <param name="offconstant">Off-center constant.</param>
        /// <returns>Coordinates of the circumcenter (or off-center)</returns>
        public Point FindCircumcenter(Point org, Point dest, Point apex,
            ref float xi, ref float eta, float offconstant)
        {
            float xdo, ydo, xao, yao;
            float dodist, aodist, dadist;
            float denominator;
            float dx, dy, dxoff, dyoff;

            Statistic.CircumcenterCount++;

            // Compute the circumcenter of the triangle.
            xdo = dest.x - org.x;
            ydo = dest.y - org.y;
            xao = apex.x - org.x;
            yao = apex.y - org.y;
            dodist = xdo * xdo + ydo * ydo;
            aodist = xao * xao + yao * yao;
            dadist = (dest.x - apex.x) * (dest.x - apex.x) +
                     (dest.y - apex.y) * (dest.y - apex.y);

            if (Behavior.NoExact)
            {
                denominator = 0.5f / (xdo * yao - xao * ydo);
            }
            else
            {
                // Use the counterclockwise() routine to ensure a positive (and
                // reasonably accurate) result, avoiding any possibility of
                // division by zero.
                denominator = 0.5f / CounterClockwise(dest, apex, org);
                // Don't count the above as an orientation test.
                Statistic.CounterClockwiseCount--;
            }

            dx = (yao * dodist - ydo * aodist) * denominator;
            dy = (xdo * aodist - xao * dodist) * denominator;

            // Find the (squared) length of the triangle's shortest edge.  This
            // serves as a conservative estimate of the insertion radius of the
            // circumcenter's parent. The estimate is used to ensure that
            // the algorithm terminates even if very small angles appear in
            // the input PSLG.
            if ((dodist < aodist) && (dodist < dadist))
            {
                if (offconstant > 0.0f)
                {
                    // Find the position of the off-center, as described by Alper Ungor.
                    dxoff = 0.5f * xdo - offconstant * ydo;
                    dyoff = 0.5f * ydo + offconstant * xdo;
                    // If the off-center is closer to the origin than the
                    // circumcenter, use the off-center instead.
                    if (dxoff * dxoff + dyoff * dyoff < dx * dx + dy * dy)
                    {
                        dx = dxoff;
                        dy = dyoff;
                    }
                }
            }
            else if (aodist < dadist)
            {
                if (offconstant > 0.0f)
                {
                    dxoff = 0.5f * xao + offconstant * yao;
                    dyoff = 0.5f * yao - offconstant * xao;
                    // If the off-center is closer to the origin than the
                    // circumcenter, use the off-center instead.
                    if (dxoff * dxoff + dyoff * dyoff < dx * dx + dy * dy)
                    {
                        dx = dxoff;
                        dy = dyoff;
                    }
                }
            }
            else
            {
                if (offconstant > 0.0f)
                {
                    dxoff = 0.5f * (apex.x - dest.x) - offconstant * (apex.y - dest.y);
                    dyoff = 0.5f * (apex.y - dest.y) + offconstant * (apex.x - dest.x);
                    // If the off-center is closer to the destination than the
                    // circumcenter, use the off-center instead.
                    if (dxoff * dxoff + dyoff * dyoff <
                        (dx - xdo) * (dx - xdo) + (dy - ydo) * (dy - ydo))
                    {
                        dx = xdo + dxoff;
                        dy = ydo + dyoff;
                    }
                }
            }

            // To interpolate vertex attributes for the new vertex inserted at
            // the circumcenter, define a coordinate system with a xi-axis,
            // directed from the triangle's origin to its destination, and
            // an eta-axis, directed from its origin to its apex.
            // Calculate the xi and eta coordinates of the circumcenter.
            xi = (yao * dx - xao * dy) * (2.0f * denominator);
            eta = (xdo * dy - ydo * dx) * (2.0f * denominator);

            return new Point(org.x + dx, org.y + dy);
        }

        /// <summary>
        /// Find the circumcenter of a triangle.
        /// </summary>
        /// <param name="org">Triangle point.</param>
        /// <param name="dest">Triangle point.</param>
        /// <param name="apex">Triangle point.</param>
        /// <param name="xi">Relative coordinate of new location.</param>
        /// <param name="eta">Relative coordinate of new location.</param>
        /// <returns>Coordinates of the circumcenter</returns>
        /// <remarks>
        /// The result is returned both in terms of x-y coordinates and xi-eta
        /// (barycentric) coordinates. The xi-eta coordinate system is defined in
        /// terms of the triangle: the origin of the triangle is the origin of the
        /// coordinate system; the destination of the triangle is one unit along the
        /// xi axis; and the apex of the triangle is one unit along the eta axis.
        /// This procedure also returns the square of the length of the triangle's
        /// shortest edge.
        /// </remarks>
        public Point FindCircumcenter(Point org, Point dest, Point apex,
            ref float xi, ref float eta)
        {
            float xdo, ydo, xao, yao;
            float dodist, aodist;
            float denominator;
            float dx, dy;

            Statistic.CircumcenterCount++;

            // Compute the circumcenter of the triangle.
            xdo = dest.x - org.x;
            ydo = dest.y - org.y;
            xao = apex.x - org.x;
            yao = apex.y - org.y;
            dodist = xdo * xdo + ydo * ydo;
            aodist = xao * xao + yao * yao;

            if (Behavior.NoExact)
            {
                denominator = 0.5f / (xdo * yao - xao * ydo);
            }
            else
            {
                // Use the counterclockwise() routine to ensure a positive (and
                // reasonably accurate) result, avoiding any possibility of
                // division by zero.
                denominator = 0.5f / CounterClockwise(dest, apex, org);
                // Don't count the above as an orientation test.
                Statistic.CounterClockwiseCount--;
            }

            dx = (yao * dodist - ydo * aodist) * denominator;
            dy = (xdo * aodist - xao * dodist) * denominator;

            // To interpolate vertex attributes for the new vertex inserted at
            // the circumcenter, define a coordinate system with a xi-axis,
            // directed from the triangle's origin to its destination, and
            // an eta-axis, directed from its origin to its apex.
            // Calculate the xi and eta coordinates of the circumcenter.
            xi = (yao * dx - xao * dy) * (2.0f * denominator);
            eta = (xdo * dy - ydo * dx) * (2.0f * denominator);

            return new Point(org.x + dx, org.y + dy);
        }

        #region Exact arithmetics

        /// <summary>
        /// Sum two expansions, eliminating zero components from the output expansion.  
        /// </summary>
        /// <param name="elen"></param>
        /// <param name="e"></param>
        /// <param name="flen"></param>
        /// <param name="f"></param>
        /// <param name="h"></param>
        /// <returns></returns>
        /// <remarks>
        /// Sets h = e + f.  See the Robust Predicates paper for details.
        /// 
        /// If round-to-even is used (as with IEEE 754), maintains the strongly nonoverlapping
        /// property.  (That is, if e is strongly nonoverlapping, h will be also.) Does NOT
        /// maintain the nonoverlapping or nonadjacent properties. 
        /// </remarks>
        private int FastExpansionSumZeroElim(int elen, float[] e, int flen, float[] f, float[] h)
        {
            float Q;
            float Qnew;
            float hh;
            float bvirt;
            float avirt, bround, around;
            int eindex, findex, hindex;
            float enow, fnow;

            enow = e[0];
            fnow = f[0];
            eindex = findex = 0;
            if ((fnow > enow) == (fnow > -enow))
            {
                Q = enow;
                enow = e[++eindex];
            }
            else
            {
                Q = fnow;
                fnow = f[++findex];
            }
            hindex = 0;
            if ((eindex < elen) && (findex < flen))
            {
                if ((fnow > enow) == (fnow > -enow))
                {
                    Qnew = (float)(enow + Q); bvirt = Qnew - enow; hh = Q - bvirt;
                    enow = e[++eindex];
                }
                else
                {
                    Qnew = (float)(fnow + Q); bvirt = Qnew - fnow; hh = Q - bvirt;
                    fnow = f[++findex];
                }
                Q = Qnew;
                if (hh != 0.0f)
                {
                    h[hindex++] = hh;
                }
                while ((eindex < elen) && (findex < flen))
                {
                    if ((fnow > enow) == (fnow > -enow))
                    {
                        Qnew = (float)(Q + enow);
                        bvirt = (float)(Qnew - Q);
                        avirt = Qnew - bvirt;
                        bround = enow - bvirt;
                        around = Q - avirt;
                        hh = around + bround;

                        enow = e[++eindex];
                    }
                    else
                    {
                        Qnew = (float)(Q + fnow);
                        bvirt = (float)(Qnew - Q);
                        avirt = Qnew - bvirt;
                        bround = fnow - bvirt;
                        around = Q - avirt;
                        hh = around + bround;

                        fnow = f[++findex];
                    }
                    Q = Qnew;
                    if (hh != 0.0f)
                    {
                        h[hindex++] = hh;
                    }
                }
            }
            while (eindex < elen)
            {
                Qnew = (float)(Q + enow);
                bvirt = (float)(Qnew - Q);
                avirt = Qnew - bvirt;
                bround = enow - bvirt;
                around = Q - avirt;
                hh = around + bround;

                enow = e[++eindex];
                Q = Qnew;
                if (hh != 0.0)
                {
                    h[hindex++] = hh;
                }
            }
            while (findex < flen)
            {
                Qnew = (float)(Q + fnow);
                bvirt = (float)(Qnew - Q);
                avirt = Qnew - bvirt;
                bround = fnow - bvirt;
                around = Q - avirt;
                hh = around + bround;

                fnow = f[++findex];
                Q = Qnew;
                if (hh != 0.0f)
                {
                    h[hindex++] = hh;
                }
            }
            if ((Q != 0.0f) || (hindex == 0))
            {
                h[hindex++] = Q;
            }
            return hindex;
        }

        /// <summary>
        /// Multiply an expansion by a scalar, eliminating zero components from the output expansion.  
        /// </summary>
        /// <param name="elen"></param>
        /// <param name="e"></param>
        /// <param name="b"></param>
        /// <param name="h"></param>
        /// <returns></returns>
        /// <remarks>
        /// Sets h = be.  See my Robust Predicates paper for details.
        /// 
        /// Maintains the nonoverlapping property.  If round-to-even is used (as with IEEE 754),
        /// maintains the strongly nonoverlapping and nonadjacent properties as well. (That is,
        /// if e has one of these properties, so will h.)
        /// </remarks>
        private int ScaleExpansionZeroElim(int elen, float[] e, float b, float[] h)
        {
            float Q, sum;
            float hh;
            float product1;
            float product0;
            int eindex, hindex;
            float enow;
            float bvirt;
            float avirt, bround, around;
            float c;
            float abig;
            float ahi, alo, bhi, blo;
            float err1, err2, err3;

            c = (float)(splitter * b); abig = (float)(c - b); bhi = c - abig; blo = b - bhi;
            Q = (float)(e[0] * b); c = (float)(splitter * e[0]); abig = (float)(c - e[0]); ahi = c - abig; alo = e[0] - ahi; err1 = Q - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); hh = (alo * blo) - err3;
            hindex = 0;
            if (hh != 0)
            {
                h[hindex++] = hh;
            }
            for (eindex = 1; eindex < elen; eindex++)
            {
                enow = e[eindex];
                product1 = (float)(enow * b); c = (float)(splitter * enow); abig = (float)(c - enow); ahi = c - abig; alo = enow - ahi; err1 = product1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); product0 = (alo * blo) - err3;
                sum = (float)(Q + product0); bvirt = (float)(sum - Q); avirt = sum - bvirt; bround = product0 - bvirt; around = Q - avirt; hh = around + bround;
                if (hh != 0)
                {
                    h[hindex++] = hh;
                }
                Q = (float)(product1 + sum); bvirt = Q - product1; hh = sum - bvirt;
                if (hh != 0)
                {
                    h[hindex++] = hh;
                }
            }
            if ((Q != 0.0f) || (hindex == 0))
            {
                h[hindex++] = Q;
            }
            return hindex;
        }

        /// <summary>
        /// Produce a one-word estimate of an expansion's value. 
        /// </summary>
        /// <param name="elen"></param>
        /// <param name="e"></param>
        /// <returns></returns>
        private float Estimate(int elen, float[] e)
        {
            float Q;
            int eindex;

            Q = e[0];
            for (eindex = 1; eindex < elen; eindex++)
            {
                Q += e[eindex];
            }
            return Q;
        }

        /// <summary>
        /// Return a positive value if the points pa, pb, and pc occur in counterclockwise
        /// order; a negative value if they occur in clockwise order; and zero if they are
        /// collinear. The result is also a rough approximation of twice the signed area of
        /// the triangle defined by the three points. 
        /// </summary>
        /// <param name="pa"></param>
        /// <param name="pb"></param>
        /// <param name="pc"></param>
        /// <param name="detsum"></param>
        /// <returns></returns>
        /// <remarks>
        /// Uses exact arithmetic if necessary to ensure a correct answer. The result returned
        /// is the determinant of a matrix. This determinant is computed adaptively, in the
        /// sense that exact arithmetic is used only to the degree it is needed to ensure that
        /// the returned value has the correct sign.  Hence, this function is usually quite fast,
        /// but will run more slowly when the input points are collinear or nearly so.
        /// </remarks>
        private float CounterClockwiseAdapt(Point pa, Point pb, Point pc, float detsum)
        {
            float acx, acy, bcx, bcy;
            float acxtail, acytail, bcxtail, bcytail;
            float detleft, detright;
            float detlefttail, detrighttail;
            float det, errbound;
            // Edited to work around index out of range exceptions (changed array length from 4 to 5).
            // See unsafe indexing in FastExpansionSumZeroElim.
            float[] B = new float[5], u = new float[5];
            float[] C1 = new float[8], C2 = new float[12], D = new float[16];
            float B3;
            int C1length, C2length, Dlength;

            float u3;
            float s1, t1;
            float s0, t0;

            float bvirt;
            float avirt, bround, around;
            float c;
            float abig;
            float ahi, alo, bhi, blo;
            float err1, err2, err3;
            float _i, _j;
            float _0;

            acx = (float)(pa.x - pc.x);
            bcx = (float)(pb.x - pc.x);
            acy = (float)(pa.y - pc.y);
            bcy = (float)(pb.y - pc.y);

            detleft = (float)(acx * bcy); c = (float)(splitter * acx); abig = (float)(c - acx); ahi = c - abig; alo = acx - ahi; c = (float)(splitter * bcy); abig = (float)(c - bcy); bhi = c - abig; blo = bcy - bhi; err1 = detleft - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); detlefttail = (alo * blo) - err3;
            detright = (float)(acy * bcx); c = (float)(splitter * acy); abig = (float)(c - acy); ahi = c - abig; alo = acy - ahi; c = (float)(splitter * bcx); abig = (float)(c - bcx); bhi = c - abig; blo = bcx - bhi; err1 = detright - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); detrighttail = (alo * blo) - err3;

            _i = (float)(detlefttail - detrighttail); bvirt = (float)(detlefttail - _i); avirt = _i + bvirt; bround = bvirt - detrighttail; around = detlefttail - avirt; B[0] = around + bround; _j = (float)(detleft + _i); bvirt = (float)(_j - detleft); avirt = _j - bvirt; bround = _i - bvirt; around = detleft - avirt; _0 = around + bround; _i = (float)(_0 - detright); bvirt = (float)(_0 - _i); avirt = _i + bvirt; bround = bvirt - detright; around = _0 - avirt; B[1] = around + bround; B3 = (float)(_j + _i); bvirt = (float)(B3 - _j); avirt = B3 - bvirt; bround = _i - bvirt; around = _j - avirt; B[2] = around + bround;

            B[3] = B3;

            det = Estimate(4, B);
            errbound = ccwerrboundB * detsum;
            if ((det >= errbound) || (-det >= errbound))
            {
                return det;
            }

            bvirt = (float)(pa.x - acx); avirt = acx + bvirt; bround = bvirt - pc.x; around = pa.x - avirt; acxtail = around + bround;
            bvirt = (float)(pb.x - bcx); avirt = bcx + bvirt; bround = bvirt - pc.x; around = pb.x - avirt; bcxtail = around + bround;
            bvirt = (float)(pa.y - acy); avirt = acy + bvirt; bround = bvirt - pc.y; around = pa.y - avirt; acytail = around + bround;
            bvirt = (float)(pb.y - bcy); avirt = bcy + bvirt; bround = bvirt - pc.y; around = pb.y - avirt; bcytail = around + bround;

            if ((acxtail == 0.0f) && (acytail == 0.0f)
                && (bcxtail == 0.0f) && (bcytail == 0.0f))
            {
                return det;
            }

            errbound = ccwerrboundC * detsum + resulterrbound * ((det) >= 0.0f ? (det) : -(det));
            det += (acx * bcytail + bcy * acxtail)
                 - (acy * bcxtail + bcx * acytail);
            if ((det >= errbound) || (-det >= errbound))
            {
                return det;
            }

            s1 = (float)(acxtail * bcy); c = (float)(splitter * acxtail); abig = (float)(c - acxtail); ahi = c - abig; alo = acxtail - ahi; c = (float)(splitter * bcy); abig = (float)(c - bcy); bhi = c - abig; blo = bcy - bhi; err1 = s1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); s0 = (alo * blo) - err3;
            t1 = (float)(acytail * bcx); c = (float)(splitter * acytail); abig = (float)(c - acytail); ahi = c - abig; alo = acytail - ahi; c = (float)(splitter * bcx); abig = (float)(c - bcx); bhi = c - abig; blo = bcx - bhi; err1 = t1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); t0 = (alo * blo) - err3;
            _i = (float)(s0 - t0); bvirt = (float)(s0 - _i); avirt = _i + bvirt; bround = bvirt - t0; around = s0 - avirt; u[0] = around + bround; _j = (float)(s1 + _i); bvirt = (float)(_j - s1); avirt = _j - bvirt; bround = _i - bvirt; around = s1 - avirt; _0 = around + bround; _i = (float)(_0 - t1); bvirt = (float)(_0 - _i); avirt = _i + bvirt; bround = bvirt - t1; around = _0 - avirt; u[1] = around + bround; u3 = (float)(_j + _i); bvirt = (float)(u3 - _j); avirt = u3 - bvirt; bround = _i - bvirt; around = _j - avirt; u[2] = around + bround;
            u[3] = u3;
            C1length = FastExpansionSumZeroElim(4, B, 4, u, C1);

            s1 = (float)(acx * bcytail); c = (float)(splitter * acx); abig = (float)(c - acx); ahi = c - abig; alo = acx - ahi; c = (float)(splitter * bcytail); abig = (float)(c - bcytail); bhi = c - abig; blo = bcytail - bhi; err1 = s1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); s0 = (alo * blo) - err3;
            t1 = (float)(acy * bcxtail); c = (float)(splitter * acy); abig = (float)(c - acy); ahi = c - abig; alo = acy - ahi; c = (float)(splitter * bcxtail); abig = (float)(c - bcxtail); bhi = c - abig; blo = bcxtail - bhi; err1 = t1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); t0 = (alo * blo) - err3;
            _i = (float)(s0 - t0); bvirt = (float)(s0 - _i); avirt = _i + bvirt; bround = bvirt - t0; around = s0 - avirt; u[0] = around + bround; _j = (float)(s1 + _i); bvirt = (float)(_j - s1); avirt = _j - bvirt; bround = _i - bvirt; around = s1 - avirt; _0 = around + bround; _i = (float)(_0 - t1); bvirt = (float)(_0 - _i); avirt = _i + bvirt; bround = bvirt - t1; around = _0 - avirt; u[1] = around + bround; u3 = (float)(_j + _i); bvirt = (float)(u3 - _j); avirt = u3 - bvirt; bround = _i - bvirt; around = _j - avirt; u[2] = around + bround;
            u[3] = u3;
            C2length = FastExpansionSumZeroElim(C1length, C1, 4, u, C2);

            s1 = (float)(acxtail * bcytail); c = (float)(splitter * acxtail); abig = (float)(c - acxtail); ahi = c - abig; alo = acxtail - ahi; c = (float)(splitter * bcytail); abig = (float)(c - bcytail); bhi = c - abig; blo = bcytail - bhi; err1 = s1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); s0 = (alo * blo) - err3;
            t1 = (float)(acytail * bcxtail); c = (float)(splitter * acytail); abig = (float)(c - acytail); ahi = c - abig; alo = acytail - ahi; c = (float)(splitter * bcxtail); abig = (float)(c - bcxtail); bhi = c - abig; blo = bcxtail - bhi; err1 = t1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); t0 = (alo * blo) - err3;
            _i = (float)(s0 - t0); bvirt = (float)(s0 - _i); avirt = _i + bvirt; bround = bvirt - t0; around = s0 - avirt; u[0] = around + bround; _j = (float)(s1 + _i); bvirt = (float)(_j - s1); avirt = _j - bvirt; bround = _i - bvirt; around = s1 - avirt; _0 = around + bround; _i = (float)(_0 - t1); bvirt = (float)(_0 - _i); avirt = _i + bvirt; bround = bvirt - t1; around = _0 - avirt; u[1] = around + bround; u3 = (float)(_j + _i); bvirt = (float)(u3 - _j); avirt = u3 - bvirt; bround = _i - bvirt; around = _j - avirt; u[2] = around + bround;
            u[3] = u3;
            Dlength = FastExpansionSumZeroElim(C2length, C2, 4, u, D);

            return (D[Dlength - 1]);
        }

        /// <summary>
        /// Return a positive value if the point pd lies inside the circle passing through
        /// pa, pb, and pc; a negative value if it lies outside; and zero if the four points
        /// are cocircular. The points pa, pb, and pc must be in counterclockwise order, or 
        /// the sign of the result will be reversed.
        /// </summary>
        /// <param name="pa"></param>
        /// <param name="pb"></param>
        /// <param name="pc"></param>
        /// <param name="pd"></param>
        /// <param name="permanent"></param>
        /// <returns></returns>
        /// <remarks>
        /// Uses exact arithmetic if necessary to ensure a correct answer. The result returned
        /// is the determinant of a matrix. This determinant is computed adaptively, in the
        /// sense that exact arithmetic is used only to the degree it is needed to ensure that
        /// the returned value has the correct sign. Hence, this function is usually quite fast,
        /// but will run more slowly when the input points are cocircular or nearly so.
        /// </remarks>
        private float InCircleAdapt(Point pa, Point pb, Point pc, Point pd, float permanent)
        {
            float adx, bdx, cdx, ady, bdy, cdy;
            float det, errbound;

            float bdxcdy1, cdxbdy1, cdxady1, adxcdy1, adxbdy1, bdxady1;
            float bdxcdy0, cdxbdy0, cdxady0, adxcdy0, adxbdy0, bdxady0;
            float[] bc = new float[4], ca = new float[4], ab = new float[4];
            float bc3, ca3, ab3;
            int axbclen, axxbclen, aybclen, ayybclen, alen;
            int bxcalen, bxxcalen, bycalen, byycalen, blen;
            int cxablen, cxxablen, cyablen, cyyablen, clen;
            int ablen;
            float[] finnow, finother, finswap;
            int finlength;

            float adxtail, bdxtail, cdxtail, adytail, bdytail, cdytail;
            float adxadx1, adyady1, bdxbdx1, bdybdy1, cdxcdx1, cdycdy1;
            float adxadx0, adyady0, bdxbdx0, bdybdy0, cdxcdx0, cdycdy0;
            float[] aa = new float[4], bb = new float[4], cc = new float[4];
            float aa3, bb3, cc3;
            float ti1, tj1;
            float ti0, tj0;
            // Edited to work around index out of range exceptions (changed array length from 4 to 5).
            // See unsafe indexing in FastExpansionSumZeroElim.
            float[] u = new float[5], v = new float[5];
            float u3, v3;
            int temp8len, temp16alen, temp16blen, temp16clen;
            int temp32alen, temp32blen, temp48len, temp64len;
            float[] axtbb = new float[8], axtcc = new float[8], aytbb = new float[8], aytcc = new float[8];
            int axtbblen, axtcclen, aytbblen, aytcclen;
            float[] bxtaa = new float[8], bxtcc = new float[8], bytaa = new float[8], bytcc = new float[8];
            int bxtaalen, bxtcclen, bytaalen, bytcclen;
            float[] cxtaa = new float[8], cxtbb = new float[8], cytaa = new float[8], cytbb = new float[8];
            int cxtaalen, cxtbblen, cytaalen, cytbblen;
            float[] axtbc = new float[8], aytbc = new float[8], bxtca = new float[8], bytca = new float[8], cxtab = new float[8], cytab = new float[8];
            int axtbclen = 0, aytbclen = 0, bxtcalen = 0, bytcalen = 0, cxtablen = 0, cytablen = 0;
            float[] axtbct = new float[16], aytbct = new float[16], bxtcat = new float[16], bytcat = new float[16], cxtabt = new float[16], cytabt = new float[16];
            int axtbctlen, aytbctlen, bxtcatlen, bytcatlen, cxtabtlen, cytabtlen;
            float[] axtbctt = new float[8], aytbctt = new float[8], bxtcatt = new float[8];
            float[] bytcatt = new float[8], cxtabtt = new float[8], cytabtt = new float[8];
            int axtbcttlen, aytbcttlen, bxtcattlen, bytcattlen, cxtabttlen, cytabttlen;
            float[] abt = new float[8], bct = new float[8], cat = new float[8];
            int abtlen, bctlen, catlen;
            float[] abtt = new float[4], bctt = new float[4], catt = new float[4];
            int abttlen, bcttlen, cattlen;
            float abtt3, bctt3, catt3;
            float negate;

            float bvirt;
            float avirt, bround, around;
            float c;
            float abig;
            float ahi, alo, bhi, blo;
            float err1, err2, err3;
            float _i, _j;
            float _0;

            adx = (float)(pa.x - pd.x);
            bdx = (float)(pb.x - pd.x);
            cdx = (float)(pc.x - pd.x);
            ady = (float)(pa.y - pd.y);
            bdy = (float)(pb.y - pd.y);
            cdy = (float)(pc.y - pd.y);

            adx = (float)(pa.x - pd.x);
            bdx = (float)(pb.x - pd.x);
            cdx = (float)(pc.x - pd.x);
            ady = (float)(pa.y - pd.y);
            bdy = (float)(pb.y - pd.y);
            cdy = (float)(pc.y - pd.y);

            bdxcdy1 = (float)(bdx * cdy); c = (float)(splitter * bdx); abig = (float)(c - bdx); ahi = c - abig; alo = bdx - ahi; c = (float)(splitter * cdy); abig = (float)(c - cdy); bhi = c - abig; blo = cdy - bhi; err1 = bdxcdy1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); bdxcdy0 = (alo * blo) - err3;
            cdxbdy1 = (float)(cdx * bdy); c = (float)(splitter * cdx); abig = (float)(c - cdx); ahi = c - abig; alo = cdx - ahi; c = (float)(splitter * bdy); abig = (float)(c - bdy); bhi = c - abig; blo = bdy - bhi; err1 = cdxbdy1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); cdxbdy0 = (alo * blo) - err3;
            _i = (float)(bdxcdy0 - cdxbdy0); bvirt = (float)(bdxcdy0 - _i); avirt = _i + bvirt; bround = bvirt - cdxbdy0; around = bdxcdy0 - avirt; bc[0] = around + bround; _j = (float)(bdxcdy1 + _i); bvirt = (float)(_j - bdxcdy1); avirt = _j - bvirt; bround = _i - bvirt; around = bdxcdy1 - avirt; _0 = around + bround; _i = (float)(_0 - cdxbdy1); bvirt = (float)(_0 - _i); avirt = _i + bvirt; bround = bvirt - cdxbdy1; around = _0 - avirt; bc[1] = around + bround; bc3 = (float)(_j + _i); bvirt = (float)(bc3 - _j); avirt = bc3 - bvirt; bround = _i - bvirt; around = _j - avirt; bc[2] = around + bround;
            bc[3] = bc3;
            axbclen = ScaleExpansionZeroElim(4, bc, adx, axbc);
            axxbclen = ScaleExpansionZeroElim(axbclen, axbc, adx, axxbc);
            aybclen = ScaleExpansionZeroElim(4, bc, ady, aybc);
            ayybclen = ScaleExpansionZeroElim(aybclen, aybc, ady, ayybc);
            alen = FastExpansionSumZeroElim(axxbclen, axxbc, ayybclen, ayybc, adet);

            cdxady1 = (float)(cdx * ady); c = (float)(splitter * cdx); abig = (float)(c - cdx); ahi = c - abig; alo = cdx - ahi; c = (float)(splitter * ady); abig = (float)(c - ady); bhi = c - abig; blo = ady - bhi; err1 = cdxady1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); cdxady0 = (alo * blo) - err3;
            adxcdy1 = (float)(adx * cdy); c = (float)(splitter * adx); abig = (float)(c - adx); ahi = c - abig; alo = adx - ahi; c = (float)(splitter * cdy); abig = (float)(c - cdy); bhi = c - abig; blo = cdy - bhi; err1 = adxcdy1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); adxcdy0 = (alo * blo) - err3;
            _i = (float)(cdxady0 - adxcdy0); bvirt = (float)(cdxady0 - _i); avirt = _i + bvirt; bround = bvirt - adxcdy0; around = cdxady0 - avirt; ca[0] = around + bround; _j = (float)(cdxady1 + _i); bvirt = (float)(_j - cdxady1); avirt = _j - bvirt; bround = _i - bvirt; around = cdxady1 - avirt; _0 = around + bround; _i = (float)(_0 - adxcdy1); bvirt = (float)(_0 - _i); avirt = _i + bvirt; bround = bvirt - adxcdy1; around = _0 - avirt; ca[1] = around + bround; ca3 = (float)(_j + _i); bvirt = (float)(ca3 - _j); avirt = ca3 - bvirt; bround = _i - bvirt; around = _j - avirt; ca[2] = around + bround;
            ca[3] = ca3;
            bxcalen = ScaleExpansionZeroElim(4, ca, bdx, bxca);
            bxxcalen = ScaleExpansionZeroElim(bxcalen, bxca, bdx, bxxca);
            bycalen = ScaleExpansionZeroElim(4, ca, bdy, byca);
            byycalen = ScaleExpansionZeroElim(bycalen, byca, bdy, byyca);
            blen = FastExpansionSumZeroElim(bxxcalen, bxxca, byycalen, byyca, bdet);

            adxbdy1 = (float)(adx * bdy); c = (float)(splitter * adx); abig = (float)(c - adx); ahi = c - abig; alo = adx - ahi; c = (float)(splitter * bdy); abig = (float)(c - bdy); bhi = c - abig; blo = bdy - bhi; err1 = adxbdy1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); adxbdy0 = (alo * blo) - err3;
            bdxady1 = (float)(bdx * ady); c = (float)(splitter * bdx); abig = (float)(c - bdx); ahi = c - abig; alo = bdx - ahi; c = (float)(splitter * ady); abig = (float)(c - ady); bhi = c - abig; blo = ady - bhi; err1 = bdxady1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); bdxady0 = (alo * blo) - err3;
            _i = (float)(adxbdy0 - bdxady0); bvirt = (float)(adxbdy0 - _i); avirt = _i + bvirt; bround = bvirt - bdxady0; around = adxbdy0 - avirt; ab[0] = around + bround; _j = (float)(adxbdy1 + _i); bvirt = (float)(_j - adxbdy1); avirt = _j - bvirt; bround = _i - bvirt; around = adxbdy1 - avirt; _0 = around + bround; _i = (float)(_0 - bdxady1); bvirt = (float)(_0 - _i); avirt = _i + bvirt; bround = bvirt - bdxady1; around = _0 - avirt; ab[1] = around + bround; ab3 = (float)(_j + _i); bvirt = (float)(ab3 - _j); avirt = ab3 - bvirt; bround = _i - bvirt; around = _j - avirt; ab[2] = around + bround;
            ab[3] = ab3;
            cxablen = ScaleExpansionZeroElim(4, ab, cdx, cxab);
            cxxablen = ScaleExpansionZeroElim(cxablen, cxab, cdx, cxxab);
            cyablen = ScaleExpansionZeroElim(4, ab, cdy, cyab);
            cyyablen = ScaleExpansionZeroElim(cyablen, cyab, cdy, cyyab);
            clen = FastExpansionSumZeroElim(cxxablen, cxxab, cyyablen, cyyab, cdet);

            ablen = FastExpansionSumZeroElim(alen, adet, blen, bdet, abdet);
            finlength = FastExpansionSumZeroElim(ablen, abdet, clen, cdet, fin1);

            det = Estimate(finlength, fin1);
            errbound = iccerrboundB * permanent;
            if ((det >= errbound) || (-det >= errbound))
            {
                return det;
            }

            bvirt = (float)(pa.x - adx); avirt = adx + bvirt; bround = bvirt - pd.x; around = pa.x - avirt; adxtail = around + bround;
            bvirt = (float)(pa.y - ady); avirt = ady + bvirt; bround = bvirt - pd.y; around = pa.y - avirt; adytail = around + bround;
            bvirt = (float)(pb.x - bdx); avirt = bdx + bvirt; bround = bvirt - pd.x; around = pb.x - avirt; bdxtail = around + bround;
            bvirt = (float)(pb.y - bdy); avirt = bdy + bvirt; bround = bvirt - pd.y; around = pb.y - avirt; bdytail = around + bround;
            bvirt = (float)(pc.x - cdx); avirt = cdx + bvirt; bround = bvirt - pd.x; around = pc.x - avirt; cdxtail = around + bround;
            bvirt = (float)(pc.y - cdy); avirt = cdy + bvirt; bround = bvirt - pd.y; around = pc.y - avirt; cdytail = around + bround;
            if ((adxtail == 0.0f) && (bdxtail == 0.0f) && (cdxtail == 0.0f)
                && (adytail == 0.0f) && (bdytail == 0.0f) && (cdytail == 0.0f))
            {
                return det;
            }

            errbound = iccerrboundC * permanent + resulterrbound * ((det) >= 0.0f ? (det) : -(det));
            det += ((adx * adx + ady * ady) * ((bdx * cdytail + cdy * bdxtail) - (bdy * cdxtail + cdx * bdytail))
                    + 2.0f * (adx * adxtail + ady * adytail) * (bdx * cdy - bdy * cdx))
                 + ((bdx * bdx + bdy * bdy) * ((cdx * adytail + ady * cdxtail) - (cdy * adxtail + adx * cdytail))
                    + 2.0f * (bdx * bdxtail + bdy * bdytail) * (cdx * ady - cdy * adx))
                 + ((cdx * cdx + cdy * cdy) * ((adx * bdytail + bdy * adxtail) - (ady * bdxtail + bdx * adytail))
                    + 2.0f * (cdx * cdxtail + cdy * cdytail) * (adx * bdy - ady * bdx));
            if ((det >= errbound) || (-det >= errbound))
            {
                return det;
            }

            finnow = fin1;
            finother = fin2;

            if ((bdxtail != 0.0f) || (bdytail != 0.0f) || (cdxtail != 0.0f) || (cdytail != 0.0f))
            {
                adxadx1 = (float)(adx * adx); c = (float)(splitter * adx); abig = (float)(c - adx); ahi = c - abig; alo = adx - ahi; err1 = adxadx1 - (ahi * ahi); err3 = err1 - ((ahi + ahi) * alo); adxadx0 = (alo * alo) - err3;
                adyady1 = (float)(ady * ady); c = (float)(splitter * ady); abig = (float)(c - ady); ahi = c - abig; alo = ady - ahi; err1 = adyady1 - (ahi * ahi); err3 = err1 - ((ahi + ahi) * alo); adyady0 = (alo * alo) - err3;
                _i = (float)(adxadx0 + adyady0); bvirt = (float)(_i - adxadx0); avirt = _i - bvirt; bround = adyady0 - bvirt; around = adxadx0 - avirt; aa[0] = around + bround; _j = (float)(adxadx1 + _i); bvirt = (float)(_j - adxadx1); avirt = _j - bvirt; bround = _i - bvirt; around = adxadx1 - avirt; _0 = around + bround; _i = (float)(_0 + adyady1); bvirt = (float)(_i - _0); avirt = _i - bvirt; bround = adyady1 - bvirt; around = _0 - avirt; aa[1] = around + bround; aa3 = (float)(_j + _i); bvirt = (float)(aa3 - _j); avirt = aa3 - bvirt; bround = _i - bvirt; around = _j - avirt; aa[2] = around + bround;
                aa[3] = aa3;
            }
            if ((cdxtail != 0.0f) || (cdytail != 0.0f) || (adxtail != 0.0f) || (adytail != 0.0f))
            {
                bdxbdx1 = (float)(bdx * bdx); c = (float)(splitter * bdx); abig = (float)(c - bdx); ahi = c - abig; alo = bdx - ahi; err1 = bdxbdx1 - (ahi * ahi); err3 = err1 - ((ahi + ahi) * alo); bdxbdx0 = (alo * alo) - err3;
                bdybdy1 = (float)(bdy * bdy); c = (float)(splitter * bdy); abig = (float)(c - bdy); ahi = c - abig; alo = bdy - ahi; err1 = bdybdy1 - (ahi * ahi); err3 = err1 - ((ahi + ahi) * alo); bdybdy0 = (alo * alo) - err3;
                _i = (float)(bdxbdx0 + bdybdy0); bvirt = (float)(_i - bdxbdx0); avirt = _i - bvirt; bround = bdybdy0 - bvirt; around = bdxbdx0 - avirt; bb[0] = around + bround; _j = (float)(bdxbdx1 + _i); bvirt = (float)(_j - bdxbdx1); avirt = _j - bvirt; bround = _i - bvirt; around = bdxbdx1 - avirt; _0 = around + bround; _i = (float)(_0 + bdybdy1); bvirt = (float)(_i - _0); avirt = _i - bvirt; bround = bdybdy1 - bvirt; around = _0 - avirt; bb[1] = around + bround; bb3 = (float)(_j + _i); bvirt = (float)(bb3 - _j); avirt = bb3 - bvirt; bround = _i - bvirt; around = _j - avirt; bb[2] = around + bround;
                bb[3] = bb3;
            }
            if ((adxtail != 0.0f) || (adytail != 0.0f) || (bdxtail != 0.0f) || (bdytail != 0.0f))
            {
                cdxcdx1 = (float)(cdx * cdx); c = (float)(splitter * cdx); abig = (float)(c - cdx); ahi = c - abig; alo = cdx - ahi; err1 = cdxcdx1 - (ahi * ahi); err3 = err1 - ((ahi + ahi) * alo); cdxcdx0 = (alo * alo) - err3;
                cdycdy1 = (float)(cdy * cdy); c = (float)(splitter * cdy); abig = (float)(c - cdy); ahi = c - abig; alo = cdy - ahi; err1 = cdycdy1 - (ahi * ahi); err3 = err1 - ((ahi + ahi) * alo); cdycdy0 = (alo * alo) - err3;
                _i = (float)(cdxcdx0 + cdycdy0); bvirt = (float)(_i - cdxcdx0); avirt = _i - bvirt; bround = cdycdy0 - bvirt; around = cdxcdx0 - avirt; cc[0] = around + bround; _j = (float)(cdxcdx1 + _i); bvirt = (float)(_j - cdxcdx1); avirt = _j - bvirt; bround = _i - bvirt; around = cdxcdx1 - avirt; _0 = around + bround; _i = (float)(_0 + cdycdy1); bvirt = (float)(_i - _0); avirt = _i - bvirt; bround = cdycdy1 - bvirt; around = _0 - avirt; cc[1] = around + bround; cc3 = (float)(_j + _i); bvirt = (float)(cc3 - _j); avirt = cc3 - bvirt; bround = _i - bvirt; around = _j - avirt; cc[2] = around + bround;
                cc[3] = cc3;
            }

            if (adxtail != 0.0f)
            {
                axtbclen = ScaleExpansionZeroElim(4, bc, adxtail, axtbc);
                temp16alen = ScaleExpansionZeroElim(axtbclen, axtbc, 2.0f * adx, temp16a);

                axtcclen = ScaleExpansionZeroElim(4, cc, adxtail, axtcc);
                temp16blen = ScaleExpansionZeroElim(axtcclen, axtcc, bdy, temp16b);

                axtbblen = ScaleExpansionZeroElim(4, bb, adxtail, axtbb);
                temp16clen = ScaleExpansionZeroElim(axtbblen, axtbb, -cdy, temp16c);

                temp32alen = FastExpansionSumZeroElim(temp16alen, temp16a, temp16blen, temp16b, temp32a);
                temp48len = FastExpansionSumZeroElim(temp16clen, temp16c, temp32alen, temp32a, temp48);
                finlength = FastExpansionSumZeroElim(finlength, finnow, temp48len, temp48, finother);
                finswap = finnow; finnow = finother; finother = finswap;
            }
            if (adytail != 0.0f)
            {
                aytbclen = ScaleExpansionZeroElim(4, bc, adytail, aytbc);
                temp16alen = ScaleExpansionZeroElim(aytbclen, aytbc, 2.0f * ady, temp16a);

                aytbblen = ScaleExpansionZeroElim(4, bb, adytail, aytbb);
                temp16blen = ScaleExpansionZeroElim(aytbblen, aytbb, cdx, temp16b);

                aytcclen = ScaleExpansionZeroElim(4, cc, adytail, aytcc);
                temp16clen = ScaleExpansionZeroElim(aytcclen, aytcc, -bdx, temp16c);

                temp32alen = FastExpansionSumZeroElim(temp16alen, temp16a, temp16blen, temp16b, temp32a);
                temp48len = FastExpansionSumZeroElim(temp16clen, temp16c, temp32alen, temp32a, temp48);
                finlength = FastExpansionSumZeroElim(finlength, finnow, temp48len, temp48, finother);
                finswap = finnow; finnow = finother; finother = finswap;
            }
            if (bdxtail != 0.0f)
            {
                bxtcalen = ScaleExpansionZeroElim(4, ca, bdxtail, bxtca);
                temp16alen = ScaleExpansionZeroElim(bxtcalen, bxtca, 2.0f * bdx, temp16a);

                bxtaalen = ScaleExpansionZeroElim(4, aa, bdxtail, bxtaa);
                temp16blen = ScaleExpansionZeroElim(bxtaalen, bxtaa, cdy, temp16b);

                bxtcclen = ScaleExpansionZeroElim(4, cc, bdxtail, bxtcc);
                temp16clen = ScaleExpansionZeroElim(bxtcclen, bxtcc, -ady, temp16c);

                temp32alen = FastExpansionSumZeroElim(temp16alen, temp16a, temp16blen, temp16b, temp32a);
                temp48len = FastExpansionSumZeroElim(temp16clen, temp16c, temp32alen, temp32a, temp48);
                finlength = FastExpansionSumZeroElim(finlength, finnow, temp48len, temp48, finother);
                finswap = finnow; finnow = finother; finother = finswap;
            }
            if (bdytail != 0.0f)
            {
                bytcalen = ScaleExpansionZeroElim(4, ca, bdytail, bytca);
                temp16alen = ScaleExpansionZeroElim(bytcalen, bytca, 2.0f * bdy, temp16a);

                bytcclen = ScaleExpansionZeroElim(4, cc, bdytail, bytcc);
                temp16blen = ScaleExpansionZeroElim(bytcclen, bytcc, adx, temp16b);

                bytaalen = ScaleExpansionZeroElim(4, aa, bdytail, bytaa);
                temp16clen = ScaleExpansionZeroElim(bytaalen, bytaa, -cdx, temp16c);

                temp32alen = FastExpansionSumZeroElim(temp16alen, temp16a, temp16blen, temp16b, temp32a);
                temp48len = FastExpansionSumZeroElim(temp16clen, temp16c, temp32alen, temp32a, temp48);
                finlength = FastExpansionSumZeroElim(finlength, finnow, temp48len, temp48, finother);
                finswap = finnow; finnow = finother; finother = finswap;
            }
            if (cdxtail != 0.0f)
            {
                cxtablen = ScaleExpansionZeroElim(4, ab, cdxtail, cxtab);
                temp16alen = ScaleExpansionZeroElim(cxtablen, cxtab, 2.0f * cdx, temp16a);

                cxtbblen = ScaleExpansionZeroElim(4, bb, cdxtail, cxtbb);
                temp16blen = ScaleExpansionZeroElim(cxtbblen, cxtbb, ady, temp16b);

                cxtaalen = ScaleExpansionZeroElim(4, aa, cdxtail, cxtaa);
                temp16clen = ScaleExpansionZeroElim(cxtaalen, cxtaa, -bdy, temp16c);

                temp32alen = FastExpansionSumZeroElim(temp16alen, temp16a, temp16blen, temp16b, temp32a);
                temp48len = FastExpansionSumZeroElim(temp16clen, temp16c, temp32alen, temp32a, temp48);
                finlength = FastExpansionSumZeroElim(finlength, finnow, temp48len, temp48, finother);
                finswap = finnow; finnow = finother; finother = finswap;
            }
            if (cdytail != 0.0f)
            {
                cytablen = ScaleExpansionZeroElim(4, ab, cdytail, cytab);
                temp16alen = ScaleExpansionZeroElim(cytablen, cytab, 2.0f * cdy, temp16a);

                cytaalen = ScaleExpansionZeroElim(4, aa, cdytail, cytaa);
                temp16blen = ScaleExpansionZeroElim(cytaalen, cytaa, bdx, temp16b);

                cytbblen = ScaleExpansionZeroElim(4, bb, cdytail, cytbb);
                temp16clen = ScaleExpansionZeroElim(cytbblen, cytbb, -adx, temp16c);

                temp32alen = FastExpansionSumZeroElim(temp16alen, temp16a, temp16blen, temp16b, temp32a);
                temp48len = FastExpansionSumZeroElim(temp16clen, temp16c, temp32alen, temp32a, temp48);
                finlength = FastExpansionSumZeroElim(finlength, finnow, temp48len, temp48, finother);
                finswap = finnow; finnow = finother; finother = finswap;
            }

            if ((adxtail != 0.0f) || (adytail != 0.0f))
            {
                if ((bdxtail != 0.0f) || (bdytail != 0.0f)
                    || (cdxtail != 0.0f) || (cdytail != 0.0f))
                {
                    ti1 = (float)(bdxtail * cdy); c = (float)(splitter * bdxtail); abig = (float)(c - bdxtail); ahi = c - abig; alo = bdxtail - ahi; c = (float)(splitter * cdy); abig = (float)(c - cdy); bhi = c - abig; blo = cdy - bhi; err1 = ti1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); ti0 = (alo * blo) - err3;
                    tj1 = (float)(bdx * cdytail); c = (float)(splitter * bdx); abig = (float)(c - bdx); ahi = c - abig; alo = bdx - ahi; c = (float)(splitter * cdytail); abig = (float)(c - cdytail); bhi = c - abig; blo = cdytail - bhi; err1 = tj1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); tj0 = (alo * blo) - err3;
                    _i = (float)(ti0 + tj0); bvirt = (float)(_i - ti0); avirt = _i - bvirt; bround = tj0 - bvirt; around = ti0 - avirt; u[0] = around + bround; _j = (float)(ti1 + _i); bvirt = (float)(_j - ti1); avirt = _j - bvirt; bround = _i - bvirt; around = ti1 - avirt; _0 = around + bround; _i = (float)(_0 + tj1); bvirt = (float)(_i - _0); avirt = _i - bvirt; bround = tj1 - bvirt; around = _0 - avirt; u[1] = around + bround; u3 = (float)(_j + _i); bvirt = (float)(u3 - _j); avirt = u3 - bvirt; bround = _i - bvirt; around = _j - avirt; u[2] = around + bround;
                    u[3] = u3;
                    negate = -bdy;
                    ti1 = (float)(cdxtail * negate); c = (float)(splitter * cdxtail); abig = (float)(c - cdxtail); ahi = c - abig; alo = cdxtail - ahi; c = (float)(splitter * negate); abig = (float)(c - negate); bhi = c - abig; blo = negate - bhi; err1 = ti1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); ti0 = (alo * blo) - err3;
                    negate = -bdytail;
                    tj1 = (float)(cdx * negate); c = (float)(splitter * cdx); abig = (float)(c - cdx); ahi = c - abig; alo = cdx - ahi; c = (float)(splitter * negate); abig = (float)(c - negate); bhi = c - abig; blo = negate - bhi; err1 = tj1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); tj0 = (alo * blo) - err3;
                    _i = (float)(ti0 + tj0); bvirt = (float)(_i - ti0); avirt = _i - bvirt; bround = tj0 - bvirt; around = ti0 - avirt; v[0] = around + bround; _j = (float)(ti1 + _i); bvirt = (float)(_j - ti1); avirt = _j - bvirt; bround = _i - bvirt; around = ti1 - avirt; _0 = around + bround; _i = (float)(_0 + tj1); bvirt = (float)(_i - _0); avirt = _i - bvirt; bround = tj1 - bvirt; around = _0 - avirt; v[1] = around + bround; v3 = (float)(_j + _i); bvirt = (float)(v3 - _j); avirt = v3 - bvirt; bround = _i - bvirt; around = _j - avirt; v[2] = around + bround;
                    v[3] = v3;
                    bctlen = FastExpansionSumZeroElim(4, u, 4, v, bct);

                    ti1 = (float)(bdxtail * cdytail); c = (float)(splitter * bdxtail); abig = (float)(c - bdxtail); ahi = c - abig; alo = bdxtail - ahi; c = (float)(splitter * cdytail); abig = (float)(c - cdytail); bhi = c - abig; blo = cdytail - bhi; err1 = ti1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); ti0 = (alo * blo) - err3;
                    tj1 = (float)(cdxtail * bdytail); c = (float)(splitter * cdxtail); abig = (float)(c - cdxtail); ahi = c - abig; alo = cdxtail - ahi; c = (float)(splitter * bdytail); abig = (float)(c - bdytail); bhi = c - abig; blo = bdytail - bhi; err1 = tj1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); tj0 = (alo * blo) - err3;
                    _i = (float)(ti0 - tj0); bvirt = (float)(ti0 - _i); avirt = _i + bvirt; bround = bvirt - tj0; around = ti0 - avirt; bctt[0] = around + bround; _j = (float)(ti1 + _i); bvirt = (float)(_j - ti1); avirt = _j - bvirt; bround = _i - bvirt; around = ti1 - avirt; _0 = around + bround; _i = (float)(_0 - tj1); bvirt = (float)(_0 - _i); avirt = _i + bvirt; bround = bvirt - tj1; around = _0 - avirt; bctt[1] = around + bround; bctt3 = (float)(_j + _i); bvirt = (float)(bctt3 - _j); avirt = bctt3 - bvirt; bround = _i - bvirt; around = _j - avirt; bctt[2] = around + bround;
                    bctt[3] = bctt3;
                    bcttlen = 4;
                }
                else
                {
                    bct[0] = 0.0f;
                    bctlen = 1;
                    bctt[0] = 0.0f;
                    bcttlen = 1;
                }

                if (adxtail != 0.0f)
                {
                    temp16alen = ScaleExpansionZeroElim(axtbclen, axtbc, adxtail, temp16a);
                    axtbctlen = ScaleExpansionZeroElim(bctlen, bct, adxtail, axtbct);
                    temp32alen = ScaleExpansionZeroElim(axtbctlen, axtbct, 2.0f * adx, temp32a);
                    temp48len = FastExpansionSumZeroElim(temp16alen, temp16a, temp32alen, temp32a, temp48);
                    finlength = FastExpansionSumZeroElim(finlength, finnow, temp48len, temp48, finother);
                    finswap = finnow; finnow = finother; finother = finswap;
                    if (bdytail != 0.0f)
                    {
                        temp8len = ScaleExpansionZeroElim(4, cc, adxtail, temp8);
                        temp16alen = ScaleExpansionZeroElim(temp8len, temp8, bdytail, temp16a);
                        finlength = FastExpansionSumZeroElim(finlength, finnow, temp16alen, temp16a, finother);
                        finswap = finnow; finnow = finother; finother = finswap;
                    }
                    if (cdytail != 0.0f)
                    {
                        temp8len = ScaleExpansionZeroElim(4, bb, -adxtail, temp8);
                        temp16alen = ScaleExpansionZeroElim(temp8len, temp8, cdytail, temp16a);
                        finlength = FastExpansionSumZeroElim(finlength, finnow, temp16alen, temp16a, finother);
                        finswap = finnow; finnow = finother; finother = finswap;
                    }

                    temp32alen = ScaleExpansionZeroElim(axtbctlen, axtbct, adxtail, temp32a);
                    axtbcttlen = ScaleExpansionZeroElim(bcttlen, bctt, adxtail, axtbctt);
                    temp16alen = ScaleExpansionZeroElim(axtbcttlen, axtbctt, 2.0f * adx, temp16a);
                    temp16blen = ScaleExpansionZeroElim(axtbcttlen, axtbctt, adxtail, temp16b);
                    temp32blen = FastExpansionSumZeroElim(temp16alen, temp16a, temp16blen, temp16b, temp32b);
                    temp64len = FastExpansionSumZeroElim(temp32alen, temp32a, temp32blen, temp32b, temp64);
                    finlength = FastExpansionSumZeroElim(finlength, finnow, temp64len, temp64, finother);
                    finswap = finnow; finnow = finother; finother = finswap;
                }
                if (adytail != 0.0)
                {
                    temp16alen = ScaleExpansionZeroElim(aytbclen, aytbc, adytail, temp16a);
                    aytbctlen = ScaleExpansionZeroElim(bctlen, bct, adytail, aytbct);
                    temp32alen = ScaleExpansionZeroElim(aytbctlen, aytbct, 2.0f * ady, temp32a);
                    temp48len = FastExpansionSumZeroElim(temp16alen, temp16a, temp32alen, temp32a, temp48);
                    finlength = FastExpansionSumZeroElim(finlength, finnow, temp48len, temp48, finother);
                    finswap = finnow; finnow = finother; finother = finswap;


                    temp32alen = ScaleExpansionZeroElim(aytbctlen, aytbct, adytail, temp32a);
                    aytbcttlen = ScaleExpansionZeroElim(bcttlen, bctt, adytail, aytbctt);
                    temp16alen = ScaleExpansionZeroElim(aytbcttlen, aytbctt, 2.0f * ady, temp16a);
                    temp16blen = ScaleExpansionZeroElim(aytbcttlen, aytbctt, adytail, temp16b);
                    temp32blen = FastExpansionSumZeroElim(temp16alen, temp16a, temp16blen, temp16b, temp32b);
                    temp64len = FastExpansionSumZeroElim(temp32alen, temp32a, temp32blen, temp32b, temp64);
                    finlength = FastExpansionSumZeroElim(finlength, finnow, temp64len, temp64, finother);
                    finswap = finnow; finnow = finother; finother = finswap;
                }
            }
            if ((bdxtail != 0.0f) || (bdytail != 0.0f))
            {
                if ((cdxtail != 0.0f) || (cdytail != 0.0f)
                    || (adxtail != 0.0f) || (adytail != 0.0f))
                {
                    ti1 = (float)(cdxtail * ady); c = (float)(splitter * cdxtail); abig = (float)(c - cdxtail); ahi = c - abig; alo = cdxtail - ahi; c = (float)(splitter * ady); abig = (float)(c - ady); bhi = c - abig; blo = ady - bhi; err1 = ti1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); ti0 = (alo * blo) - err3;
                    tj1 = (float)(cdx * adytail); c = (float)(splitter * cdx); abig = (float)(c - cdx); ahi = c - abig; alo = cdx - ahi; c = (float)(splitter * adytail); abig = (float)(c - adytail); bhi = c - abig; blo = adytail - bhi; err1 = tj1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); tj0 = (alo * blo) - err3;
                    _i = (float)(ti0 + tj0); bvirt = (float)(_i - ti0); avirt = _i - bvirt; bround = tj0 - bvirt; around = ti0 - avirt; u[0] = around + bround; _j = (float)(ti1 + _i); bvirt = (float)(_j - ti1); avirt = _j - bvirt; bround = _i - bvirt; around = ti1 - avirt; _0 = around + bround; _i = (float)(_0 + tj1); bvirt = (float)(_i - _0); avirt = _i - bvirt; bround = tj1 - bvirt; around = _0 - avirt; u[1] = around + bround; u3 = (float)(_j + _i); bvirt = (float)(u3 - _j); avirt = u3 - bvirt; bround = _i - bvirt; around = _j - avirt; u[2] = around + bround;
                    u[3] = u3;
                    negate = -cdy;
                    ti1 = (float)(adxtail * negate); c = (float)(splitter * adxtail); abig = (float)(c - adxtail); ahi = c - abig; alo = adxtail - ahi; c = (float)(splitter * negate); abig = (float)(c - negate); bhi = c - abig; blo = negate - bhi; err1 = ti1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); ti0 = (alo * blo) - err3;
                    negate = -cdytail;
                    tj1 = (float)(adx * negate); c = (float)(splitter * adx); abig = (float)(c - adx); ahi = c - abig; alo = adx - ahi; c = (float)(splitter * negate); abig = (float)(c - negate); bhi = c - abig; blo = negate - bhi; err1 = tj1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); tj0 = (alo * blo) - err3;
                    _i = (float)(ti0 + tj0); bvirt = (float)(_i - ti0); avirt = _i - bvirt; bround = tj0 - bvirt; around = ti0 - avirt; v[0] = around + bround; _j = (float)(ti1 + _i); bvirt = (float)(_j - ti1); avirt = _j - bvirt; bround = _i - bvirt; around = ti1 - avirt; _0 = around + bround; _i = (float)(_0 + tj1); bvirt = (float)(_i - _0); avirt = _i - bvirt; bround = tj1 - bvirt; around = _0 - avirt; v[1] = around + bround; v3 = (float)(_j + _i); bvirt = (float)(v3 - _j); avirt = v3 - bvirt; bround = _i - bvirt; around = _j - avirt; v[2] = around + bround;
                    v[3] = v3;
                    catlen = FastExpansionSumZeroElim(4, u, 4, v, cat);

                    ti1 = (float)(cdxtail * adytail); c = (float)(splitter * cdxtail); abig = (float)(c - cdxtail); ahi = c - abig; alo = cdxtail - ahi; c = (float)(splitter * adytail); abig = (float)(c - adytail); bhi = c - abig; blo = adytail - bhi; err1 = ti1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); ti0 = (alo * blo) - err3;
                    tj1 = (float)(adxtail * cdytail); c = (float)(splitter * adxtail); abig = (float)(c - adxtail); ahi = c - abig; alo = adxtail - ahi; c = (float)(splitter * cdytail); abig = (float)(c - cdytail); bhi = c - abig; blo = cdytail - bhi; err1 = tj1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); tj0 = (alo * blo) - err3;
                    _i = (float)(ti0 - tj0); bvirt = (float)(ti0 - _i); avirt = _i + bvirt; bround = bvirt - tj0; around = ti0 - avirt; catt[0] = around + bround; _j = (float)(ti1 + _i); bvirt = (float)(_j - ti1); avirt = _j - bvirt; bround = _i - bvirt; around = ti1 - avirt; _0 = around + bround; _i = (float)(_0 - tj1); bvirt = (float)(_0 - _i); avirt = _i + bvirt; bround = bvirt - tj1; around = _0 - avirt; catt[1] = around + bround; catt3 = (float)(_j + _i); bvirt = (float)(catt3 - _j); avirt = catt3 - bvirt; bround = _i - bvirt; around = _j - avirt; catt[2] = around + bround;
                    catt[3] = catt3;
                    cattlen = 4;
                }
                else
                {
                    cat[0] = 0.0f;
                    catlen = 1;
                    catt[0] = 0.0f;
                    cattlen = 1;
                }

                if (bdxtail != 0.0f)
                {
                    temp16alen = ScaleExpansionZeroElim(bxtcalen, bxtca, bdxtail, temp16a);
                    bxtcatlen = ScaleExpansionZeroElim(catlen, cat, bdxtail, bxtcat);
                    temp32alen = ScaleExpansionZeroElim(bxtcatlen, bxtcat, 2.0f * bdx, temp32a);
                    temp48len = FastExpansionSumZeroElim(temp16alen, temp16a, temp32alen, temp32a, temp48);
                    finlength = FastExpansionSumZeroElim(finlength, finnow, temp48len, temp48, finother);
                    finswap = finnow; finnow = finother; finother = finswap;
                    if (cdytail != 0.0f)
                    {
                        temp8len = ScaleExpansionZeroElim(4, aa, bdxtail, temp8);
                        temp16alen = ScaleExpansionZeroElim(temp8len, temp8, cdytail, temp16a);
                        finlength = FastExpansionSumZeroElim(finlength, finnow, temp16alen, temp16a, finother);
                        finswap = finnow; finnow = finother; finother = finswap;
                    }
                    if (adytail != 0.0f)
                    {
                        temp8len = ScaleExpansionZeroElim(4, cc, -bdxtail, temp8);
                        temp16alen = ScaleExpansionZeroElim(temp8len, temp8, adytail, temp16a);
                        finlength = FastExpansionSumZeroElim(finlength, finnow, temp16alen, temp16a, finother);
                        finswap = finnow; finnow = finother; finother = finswap;
                    }

                    temp32alen = ScaleExpansionZeroElim(bxtcatlen, bxtcat, bdxtail, temp32a);
                    bxtcattlen = ScaleExpansionZeroElim(cattlen, catt, bdxtail, bxtcatt);
                    temp16alen = ScaleExpansionZeroElim(bxtcattlen, bxtcatt, 2.0f * bdx, temp16a);
                    temp16blen = ScaleExpansionZeroElim(bxtcattlen, bxtcatt, bdxtail, temp16b);
                    temp32blen = FastExpansionSumZeroElim(temp16alen, temp16a, temp16blen, temp16b, temp32b);
                    temp64len = FastExpansionSumZeroElim(temp32alen, temp32a, temp32blen, temp32b, temp64);
                    finlength = FastExpansionSumZeroElim(finlength, finnow, temp64len, temp64, finother);
                    finswap = finnow; finnow = finother; finother = finswap;
                }
                if (bdytail != 0.0f)
                {
                    temp16alen = ScaleExpansionZeroElim(bytcalen, bytca, bdytail, temp16a);
                    bytcatlen = ScaleExpansionZeroElim(catlen, cat, bdytail, bytcat);
                    temp32alen = ScaleExpansionZeroElim(bytcatlen, bytcat, 2.0f * bdy, temp32a);
                    temp48len = FastExpansionSumZeroElim(temp16alen, temp16a, temp32alen, temp32a, temp48);
                    finlength = FastExpansionSumZeroElim(finlength, finnow, temp48len, temp48, finother);
                    finswap = finnow; finnow = finother; finother = finswap;

                    temp32alen = ScaleExpansionZeroElim(bytcatlen, bytcat, bdytail, temp32a);
                    bytcattlen = ScaleExpansionZeroElim(cattlen, catt, bdytail, bytcatt);
                    temp16alen = ScaleExpansionZeroElim(bytcattlen, bytcatt, 2.0f * bdy, temp16a);
                    temp16blen = ScaleExpansionZeroElim(bytcattlen, bytcatt, bdytail, temp16b);
                    temp32blen = FastExpansionSumZeroElim(temp16alen, temp16a, temp16blen, temp16b, temp32b);
                    temp64len = FastExpansionSumZeroElim(temp32alen, temp32a, temp32blen, temp32b, temp64);
                    finlength = FastExpansionSumZeroElim(finlength, finnow, temp64len, temp64, finother);
                    finswap = finnow; finnow = finother; finother = finswap;
                }
            }
            if ((cdxtail != 0.0f) || (cdytail != 0.0f))
            {
                if ((adxtail != 0.0f) || (adytail != 0.0f)
                    || (bdxtail != 0.0f) || (bdytail != 0.0f))
                {
                    ti1 = (float)(adxtail * bdy); c = (float)(splitter * adxtail); abig = (float)(c - adxtail); ahi = c - abig; alo = adxtail - ahi; c = (float)(splitter * bdy); abig = (float)(c - bdy); bhi = c - abig; blo = bdy - bhi; err1 = ti1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); ti0 = (alo * blo) - err3;
                    tj1 = (float)(adx * bdytail); c = (float)(splitter * adx); abig = (float)(c - adx); ahi = c - abig; alo = adx - ahi; c = (float)(splitter * bdytail); abig = (float)(c - bdytail); bhi = c - abig; blo = bdytail - bhi; err1 = tj1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); tj0 = (alo * blo) - err3;
                    _i = (float)(ti0 + tj0); bvirt = (float)(_i - ti0); avirt = _i - bvirt; bround = tj0 - bvirt; around = ti0 - avirt; u[0] = around + bround; _j = (float)(ti1 + _i); bvirt = (float)(_j - ti1); avirt = _j - bvirt; bround = _i - bvirt; around = ti1 - avirt; _0 = around + bround; _i = (float)(_0 + tj1); bvirt = (float)(_i - _0); avirt = _i - bvirt; bround = tj1 - bvirt; around = _0 - avirt; u[1] = around + bround; u3 = (float)(_j + _i); bvirt = (float)(u3 - _j); avirt = u3 - bvirt; bround = _i - bvirt; around = _j - avirt; u[2] = around + bround;
                    u[3] = u3;
                    negate = -ady;
                    ti1 = (float)(bdxtail * negate); c = (float)(splitter * bdxtail); abig = (float)(c - bdxtail); ahi = c - abig; alo = bdxtail - ahi; c = (float)(splitter * negate); abig = (float)(c - negate); bhi = c - abig; blo = negate - bhi; err1 = ti1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); ti0 = (alo * blo) - err3;
                    negate = -adytail;
                    tj1 = (float)(bdx * negate); c = (float)(splitter * bdx); abig = (float)(c - bdx); ahi = c - abig; alo = bdx - ahi; c = (float)(splitter * negate); abig = (float)(c - negate); bhi = c - abig; blo = negate - bhi; err1 = tj1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); tj0 = (alo * blo) - err3;
                    _i = (float)(ti0 + tj0); bvirt = (float)(_i - ti0); avirt = _i - bvirt; bround = tj0 - bvirt; around = ti0 - avirt; v[0] = around + bround; _j = (float)(ti1 + _i); bvirt = (float)(_j - ti1); avirt = _j - bvirt; bround = _i - bvirt; around = ti1 - avirt; _0 = around + bround; _i = (float)(_0 + tj1); bvirt = (float)(_i - _0); avirt = _i - bvirt; bround = tj1 - bvirt; around = _0 - avirt; v[1] = around + bround; v3 = (float)(_j + _i); bvirt = (float)(v3 - _j); avirt = v3 - bvirt; bround = _i - bvirt; around = _j - avirt; v[2] = around + bround;
                    v[3] = v3;
                    abtlen = FastExpansionSumZeroElim(4, u, 4, v, abt);

                    ti1 = (float)(adxtail * bdytail); c = (float)(splitter * adxtail); abig = (float)(c - adxtail); ahi = c - abig; alo = adxtail - ahi; c = (float)(splitter * bdytail); abig = (float)(c - bdytail); bhi = c - abig; blo = bdytail - bhi; err1 = ti1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); ti0 = (alo * blo) - err3;
                    tj1 = (float)(bdxtail * adytail); c = (float)(splitter * bdxtail); abig = (float)(c - bdxtail); ahi = c - abig; alo = bdxtail - ahi; c = (float)(splitter * adytail); abig = (float)(c - adytail); bhi = c - abig; blo = adytail - bhi; err1 = tj1 - (ahi * bhi); err2 = err1 - (alo * bhi); err3 = err2 - (ahi * blo); tj0 = (alo * blo) - err3;
                    _i = (float)(ti0 - tj0); bvirt = (float)(ti0 - _i); avirt = _i + bvirt; bround = bvirt - tj0; around = ti0 - avirt; abtt[0] = around + bround; _j = (float)(ti1 + _i); bvirt = (float)(_j - ti1); avirt = _j - bvirt; bround = _i - bvirt; around = ti1 - avirt; _0 = around + bround; _i = (float)(_0 - tj1); bvirt = (float)(_0 - _i); avirt = _i + bvirt; bround = bvirt - tj1; around = _0 - avirt; abtt[1] = around + bround; abtt3 = (float)(_j + _i); bvirt = (float)(abtt3 - _j); avirt = abtt3 - bvirt; bround = _i - bvirt; around = _j - avirt; abtt[2] = around + bround;
                    abtt[3] = abtt3;
                    abttlen = 4;
                }
                else
                {
                    abt[0] = 0.0f;
                    abtlen = 1;
                    abtt[0] = 0.0f;
                    abttlen = 1;
                }

                if (cdxtail != 0.0f)
                {
                    temp16alen = ScaleExpansionZeroElim(cxtablen, cxtab, cdxtail, temp16a);
                    cxtabtlen = ScaleExpansionZeroElim(abtlen, abt, cdxtail, cxtabt);
                    temp32alen = ScaleExpansionZeroElim(cxtabtlen, cxtabt, 2.0f * cdx, temp32a);
                    temp48len = FastExpansionSumZeroElim(temp16alen, temp16a, temp32alen, temp32a, temp48);
                    finlength = FastExpansionSumZeroElim(finlength, finnow, temp48len, temp48, finother);
                    finswap = finnow; finnow = finother; finother = finswap;
                    if (adytail != 0.0f)
                    {
                        temp8len = ScaleExpansionZeroElim(4, bb, cdxtail, temp8);
                        temp16alen = ScaleExpansionZeroElim(temp8len, temp8, adytail, temp16a);
                        finlength = FastExpansionSumZeroElim(finlength, finnow, temp16alen, temp16a, finother);
                        finswap = finnow; finnow = finother; finother = finswap;
                    }
                    if (bdytail != 0.0f)
                    {
                        temp8len = ScaleExpansionZeroElim(4, aa, -cdxtail, temp8);
                        temp16alen = ScaleExpansionZeroElim(temp8len, temp8, bdytail, temp16a);
                        finlength = FastExpansionSumZeroElim(finlength, finnow, temp16alen, temp16a, finother);
                        finswap = finnow; finnow = finother; finother = finswap;
                    }

                    temp32alen = ScaleExpansionZeroElim(cxtabtlen, cxtabt, cdxtail, temp32a);
                    cxtabttlen = ScaleExpansionZeroElim(abttlen, abtt, cdxtail, cxtabtt);
                    temp16alen = ScaleExpansionZeroElim(cxtabttlen, cxtabtt, 2.0f * cdx, temp16a);
                    temp16blen = ScaleExpansionZeroElim(cxtabttlen, cxtabtt, cdxtail, temp16b);
                    temp32blen = FastExpansionSumZeroElim(temp16alen, temp16a, temp16blen, temp16b, temp32b);
                    temp64len = FastExpansionSumZeroElim(temp32alen, temp32a, temp32blen, temp32b, temp64);
                    finlength = FastExpansionSumZeroElim(finlength, finnow, temp64len, temp64, finother);
                    finswap = finnow; finnow = finother; finother = finswap;
                }
                if (cdytail != 0.0f)
                {
                    temp16alen = ScaleExpansionZeroElim(cytablen, cytab, cdytail, temp16a);
                    cytabtlen = ScaleExpansionZeroElim(abtlen, abt, cdytail, cytabt);
                    temp32alen = ScaleExpansionZeroElim(cytabtlen, cytabt, 2.0f * cdy, temp32a);
                    temp48len = FastExpansionSumZeroElim(temp16alen, temp16a, temp32alen, temp32a, temp48);
                    finlength = FastExpansionSumZeroElim(finlength, finnow, temp48len, temp48, finother);
                    finswap = finnow; finnow = finother; finother = finswap;


                    temp32alen = ScaleExpansionZeroElim(cytabtlen, cytabt, cdytail, temp32a);
                    cytabttlen = ScaleExpansionZeroElim(abttlen, abtt, cdytail, cytabtt);
                    temp16alen = ScaleExpansionZeroElim(cytabttlen, cytabtt, 2.0f * cdy, temp16a);
                    temp16blen = ScaleExpansionZeroElim(cytabttlen, cytabtt, cdytail, temp16b);
                    temp32blen = FastExpansionSumZeroElim(temp16alen, temp16a, temp16blen, temp16b, temp32b);
                    temp64len = FastExpansionSumZeroElim(temp32alen, temp32a, temp32blen, temp32b, temp64);
                    finlength = FastExpansionSumZeroElim(finlength, finnow, temp64len, temp64, finother);
                    finswap = finnow; finnow = finother; finother = finswap;
                }
            }

            return finnow[finlength - 1];
        }

        #region Workspace

        // InCircleAdapt workspace:
        float[] fin1, fin2, abdet;

        float[] axbc, axxbc, aybc, ayybc, adet;
        float[] bxca, bxxca, byca, byyca, bdet;
        float[] cxab, cxxab, cyab, cyyab, cdet;

        float[] temp8, temp16a, temp16b, temp16c;
        float[] temp32a, temp32b, temp48, temp64;

        private void AllocateWorkspace()
        {
            fin1 = new float[1152];
            fin2 = new float[1152];
            abdet = new float[64];

            axbc = new float[8];
            axxbc = new float[16];
            aybc = new float[8];
            ayybc = new float[16];
            adet = new float[32];

            bxca = new float[8];
            bxxca = new float[16];
            byca = new float[8];
            byyca = new float[16];
            bdet = new float[32];

            cxab = new float[8];
            cxxab = new float[16];
            cyab = new float[8];
            cyyab = new float[16];
            cdet = new float[32];

            temp8 = new float[8];
            temp16a = new float[16];
            temp16b = new float[16];
            temp16c = new float[16];

            temp32a = new float[32];
            temp32b = new float[32];
            temp48 = new float[48];
            temp64 = new float[64];
        }

        private void ClearWorkspace()
        {
        }

        #endregion

        #endregion
    }
}
