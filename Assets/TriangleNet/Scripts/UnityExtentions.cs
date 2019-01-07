using System.Collections;
using System.Collections.Generic;
using TriangleNet.Geometry;
using UnityEngine;

namespace TriangleNet
{
    public static class UnityExtentions
    {
        public static Polygon FromContourToPolygon(this List<Vector2> contour, List<List<Vector2>> holes = null)
        {
            Polygon polygon = new Polygon();
  
            polygon.Add(new Contour(contour.ToTriangleNetVertices()));
            if (holes != null)
            {
                foreach (var hole in holes)
                {
                    polygon.Add(new Contour(hole.ToTriangleNetVertices()), true);
                }
            }
            
            return polygon;
        }

        public static Polygon FromRandomVertexToPolygon(this List<Vector2> vertices)
        {
            
        }
        public static List<Vertex> ToTriangleNetVertices(this List<Vector2> points)
        {
            List<Vertex> vertices = new List<Vertex>();
            foreach (var vec in points)
            {
                vertices.Add(new Vertex(vec.x, vec.y));
            }

            return vertices;
        }

        public static void DrawGizmos(this TriangleNetMesh triangleNetMesh)
        {
            foreach (var triangle in triangleNetMesh.triangles)
            {
                var verts = triangle.vertices;
                Gizmos.DrawLine((Vector3) verts[0], (Vector3) verts[1]);
                Gizmos.DrawLine((Vector3) verts[1], (Vector3) verts[2]);
                Gizmos.DrawLine((Vector3) verts[2], (Vector3) verts[0]);	
            }
        }
    }
}