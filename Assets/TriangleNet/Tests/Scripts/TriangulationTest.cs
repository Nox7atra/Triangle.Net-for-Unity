using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using TriangleNet;
using TriangleNet.Geometry;
using TriangleNet.Meshing;
using UnityEngine;
using Debug = UnityEngine.Debug;
using Random = UnityEngine.Random;

public class TriangulationTest : MonoBehaviour
{
	[SerializeField] private int _VertexCount;
	
	[SerializeField] private List<Vector2> _Contour;
	[SerializeField] private MeshFilter _Filter;
	private Polygon _Polygon;
	private TriangleNetMesh _TriangleNetMesh;

	public void GenerateRandomPoints()
	{
		_Contour.Clear();
		for (int i = 0; i < _VertexCount; i++)
		{
			_Contour.Add(new Vector2(Random.Range(0, 10f), Random.Range(0, 10f)));
		}
	}
	public void Triangulate()
	{
		_Polygon = new Polygon();
		foreach (var vector2 in _Contour)
		{
			_Polygon.Add(vector2);
		}
		_TriangleNetMesh = (TriangleNetMesh) _Polygon.Triangulate();
		_Filter.mesh = _TriangleNetMesh.GenerateUnityMesh();
	}

	private void OnDrawGizmosSelected()
	{
		if(_TriangleNetMesh != null)
			_TriangleNetMesh.DrawGizmos();
	}
}

#if UNITY_EDITOR
[UnityEditor.CustomEditor(typeof(TriangulationTest))]
public class TriangulatorTestEditor : UnityEditor.Editor
{
	public override void OnInspectorGUI()
	{
		base.OnInspectorGUI();
		if (GUILayout.Button("Generate points"))
		{
			var test = target as TriangulationTest;
			test.GenerateRandomPoints();
		}
		if (GUILayout.Button("Triangulate"))
		{
			var test = target as TriangulationTest;
			test.Triangulate();
			UnityEditor.SceneView.RepaintAll();
		}
      
	}
}
#endif