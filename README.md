# Triangle.Net-for-Unity

Установка через Unity Package Manager / Add package from git URL: https://github.com/Nox7atra/Triangle.Net-for-Unity.git?path=/Assets/TriangleNet

Article about (rus): https://habr.com/post/435374/

## Extention Methods

*public static void Add(this Polygon polygon, List<Vector2> contour, bool isHole = false)* - Add contour to Polygon

*public static void Add(this Polygon polygon, Vector2 vertex)* - Add vertex to Polygon

*public static Mesh GenerateUnityMesh(this TriangleNetMesh triangleNetMesh, QualityOptions options = null)* - Generate Unity Mesh
