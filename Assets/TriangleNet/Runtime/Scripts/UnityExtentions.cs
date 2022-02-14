using System.Collections;
using System.Collections.Generic;
using System.Linq;
using TriangleNet.Geometry;
using TriangleNet.Meshing;
using UnityEngine;

namespace TriangleNet
{
    public static class UnityExtentions
    {
        public static void Add(this Polygon polygon, List<Vector2> contour, bool isHole = false)
        {
            polygon.Add(new Contour(contour.ToTriangleNetVertices()), isHole);
        }

        public static void Add(this Polygon polygon, Vector2 vertex)
        {
            polygon.Add(new Vertex(vertex.x, vertex.y));
        }

        public static Mesh GenerateUnityMesh(this TriangleNetMesh triangleNetMesh, QualityOptions options = null)
        {
            if (options != null)
            {
                triangleNetMesh.Refine(options);
            }
         
            Mesh mesh = new Mesh();
            var triangleNetVerts = triangleNetMesh.Vertices.ToList();
  
            var triangles = triangleNetMesh.Triangles;
       
            Vector3[] verts = new Vector3[triangleNetVerts.Count];
            int[] trisIndex = new int[triangles.Count * 3];

            for (int i = 0; i < verts.Length; i++)
            {
                verts[i] = (Vector3) triangleNetVerts[i];
            }
            
            int k = 0;
         
            foreach (var triangle in triangles)
            {
                for (int i = 2; i >= 0; i--)
                {
                    trisIndex[k] = triangleNetVerts.IndexOf(triangle.GetVertex(i));
                    k++;
                }
            }

            mesh.vertices = verts;
            mesh.triangles = trisIndex;

            mesh.RecalculateBounds();
            mesh.RecalculateNormals();
            return mesh;
        }
        
        private static List<Vertex> ToTriangleNetVertices(this List<Vector2> points)
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