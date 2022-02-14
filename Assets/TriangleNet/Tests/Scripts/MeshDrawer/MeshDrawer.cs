using System.Collections;
using System.Collections.Generic;
using System.Linq;
using TriangleNet.Geometry;
using UnityEngine;

namespace TriangleNet
{
	public class MeshDrawer : MonoBehaviour
	{
		[SerializeField] private Material _MeshMaterial;
		[SerializeField] private Material _ContourMaterial;
		[SerializeField] private Material _HoleMaterial;
		[SerializeField] private float _UVScale;
		[SerializeField] private float _ContourLockRadius;
		private List<Vector2> _Contour;
		private List<List<Vector2>> _Holes;
		
		private MeshDrawerState _CurrentState;
		private List<Vector2> _CurrentHole;
		private Vector3 _CurrentMousePosition;
		
		private void Start()
		{
			_Contour = new List<Vector2>();
			_Holes = new List<List<Vector2>>();
		}

		void Update()
		{
			if(_CurrentState == MeshDrawerState.Nothing) return;;
			
			_CurrentMousePosition = Camera.main.ScreenToWorldPoint(Input.mousePosition);
			if (Input.GetMouseButtonDown(0))
			{
				var currentContour = _CurrentState == MeshDrawerState.DrawContour ? _Contour : _CurrentHole;
				if (currentContour.Count > 1 && Vector2.Distance(currentContour[0], _CurrentMousePosition) < _ContourLockRadius)
				{
					if (_CurrentState == MeshDrawerState.DrawHole)
					{
						_Holes.Add(_CurrentHole);
						_CurrentHole = null;
					}
					_CurrentState = MeshDrawerState.Nothing;
					
				}
				else
				{
					currentContour.Add(_CurrentMousePosition);
				}
			}
		}
		
		public void StartDrawContour()
		{
			if(_CurrentState != MeshDrawerState.Nothing) return;
			Clear();
			_CurrentState = MeshDrawerState.DrawContour;
		}
		public void StartDrawHole()
		{
			if(_CurrentState != MeshDrawerState.Nothing) return;
			_CurrentState = MeshDrawerState.DrawHole;
			_CurrentHole = new List<Vector2>();
		
		}
		public void GenerateMesh()
		{
			if(_CurrentState != MeshDrawerState.Nothing) return;
			Polygon poly = new Polygon();
			poly.Add(_Contour);
			foreach (var hole in _Holes)
			{
				poly.Add(hole, true);
			}
		
			var triangleNetMesh = (TriangleNetMesh) poly.Triangulate();
			
			GameObject go = new GameObject("Generated mesh");
			var mf = go.AddComponent<MeshFilter>();
			var mesh = triangleNetMesh.GenerateUnityMesh();
			mesh.uv = GenerateUv(mesh.vertices);
			mf.mesh = mesh;
			var mr = go.AddComponent<MeshRenderer>();
			mr.sharedMaterial = _MeshMaterial;
		
			var collider = go.AddComponent<PolygonCollider2D>();
			collider.points = _Contour.ToArray();
		
			var rb = go.AddComponent<Rigidbody2D>();
			rb.mass = triangleNetMesh.Triangles.Sum(tris => tris.Area);
			Clear();
		}

		private Vector2[] GenerateUv(Vector3[] vertices)
		{
			Vector2[] uvs = new Vector2[vertices.Length];
			for (int i = 0; i < vertices.Length; i++)
			{
				uvs[i] = new Vector2(vertices[i].x * _UVScale, vertices[i].y * _UVScale);
			}

			return uvs;
		}
		private void Clear()
		{
			
			_Contour.Clear();
			_Holes.Clear();
		}
		
		private void OnRenderObject()
		{
			if (_Contour.Count > 0)
			{
				_ContourMaterial.SetPass(0);
				GL.PushMatrix();
				GL.Begin(GL.LINES);
				for (int i = 0; i < _Contour.Count - 1; ++i)
				{
					GL.Vertex(_Contour[i]);
					GL.Vertex(_Contour[i + 1]);
				}

				if (_CurrentState == MeshDrawerState.DrawContour)
				{
					GL.Vertex(_Contour[_Contour.Count - 1]);
					GL.Vertex(_CurrentMousePosition);
				}
				else
				{
					GL.Vertex(_Contour[_Contour.Count - 1]);
					GL.Vertex(_Contour[0]);
				}
				GL.End();
				GL.PopMatrix();
			}
			
			_HoleMaterial.SetPass(0);
			GL.PushMatrix();
			GL.Begin(GL.LINES);			
			if (_Holes.Count > 0)
			{
				for (int i = 0; i < _Holes.Count; i++)
				{
					for (int j = 0; j < _Holes[i].Count; ++j)
					{
						GL.Vertex(_Holes[i][j]);
						GL.Vertex(_Holes[i][(j + 1) % _Holes[i].Count]);
					}
				}
			}
			if (_CurrentHole != null && _CurrentHole.Count > 0)
			{
				for (int j = 0; j < _CurrentHole.Count - 1; ++j)
				{
					GL.Vertex(_CurrentHole[j]);
					GL.Vertex(_CurrentHole[(j + 1)]);
				}
				if (_CurrentState == MeshDrawerState.DrawHole)
				{
					GL.Vertex(_CurrentHole[_CurrentHole.Count - 1]);
					GL.Vertex(_CurrentMousePosition);
				}
			}
			

			GL.End();
			GL.PopMatrix();
		}
		private enum MeshDrawerState
		{
			Nothing,
			DrawContour,
			DrawHole
		}
	}
}