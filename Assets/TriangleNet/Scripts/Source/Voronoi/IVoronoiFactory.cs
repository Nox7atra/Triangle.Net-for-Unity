
namespace TriangleNet.Voronoi
{
    using TriangleNet.Topology.DCEL;

    public interface IVoronoiFactory
    {
        void Initialize(int vertexCount, int edgeCount, int faceCount);

        void Reset();

        Vertex CreateVertex(float x, float y);

        HalfEdge CreateHalfEdge(Vertex origin, Face face);

        Face CreateFace(Geometry.Vertex vertex);
    }
}
