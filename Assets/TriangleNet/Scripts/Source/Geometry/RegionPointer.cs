// -----------------------------------------------------------------------
// <copyright file="RegionPointer.cs" company="">
// Triangle.NET code by Christian Woltering, http://triangle.codeplex.com/
// </copyright>
// -----------------------------------------------------------------------

namespace TriangleNet.Geometry
{
    using System;
    using System.Collections.Generic;

    /// <summary>
    /// Pointer to a region in the mesh geometry. A region is a well-defined
    /// subset of the geomerty (enclosed by subsegments).
    /// </summary>
    public class RegionPointer
    {
        internal Point point;
        internal int id;
        internal float area;

        /// <summary>
        /// Gets or sets a region area constraint.
        /// </summary>
        public float Area
        {
            get { return area; }
            set
            {
                if (value < 0.0f)
                {
                    throw new ArgumentException("Area constraints must not be negative.");
                }
                area = value;
            }
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="RegionPointer" /> class.
        /// </summary>
        /// <param name="x">X coordinate of the region.</param>
        /// <param name="y">Y coordinate of the region.</param>
        /// <param name="id">Region id.</param>
        public RegionPointer(float x, float y, int id)
            : this(x, y, id, 0.0f)
        {
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="RegionPointer" /> class.
        /// </summary>
        /// <param name="x">X coordinate of the region.</param>
        /// <param name="y">Y coordinate of the region.</param>
        /// <param name="id">Region id.</param>
        /// <param name="area">Area constraint.</param>
        public RegionPointer(float x, float y, int id, float area)
        {
            this.point = new Point(x, y);
            this.id = id;
            this.area = area;
        }
    }
}
