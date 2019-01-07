using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
public static class MathUtils
{
    public static bool IsPointInsideCountour(List<Vector2> contour, Vector2 point)
	{
		Vector2 p1, p2;
		bool inside = false;

		if (contour.Count < 3)
		{
			return false;
		}

		Vector2 oldPoint = new Vector2(contour[contour.Count - 1].x, contour[contour.Count - 1].y);

		for (int i = 0; i < contour.Count; i++)
		{
			Vector2 newPoint = new Vector2(contour[i].x, contour[i].y);

			if (newPoint.x > oldPoint.x)
			{
				p1 = oldPoint;
				p2 = newPoint;
			}
			else
			{
				p1 = newPoint;
				p2 = oldPoint;
			}

			if ((newPoint.x < point.x) == (point.x <= oldPoint.x) &&
			   (point.y - p1.y) * (p2.x - p1.x) < (p2.y - p1.y) * (point.x - p1.x))
			{
				inside = !inside;
			}

			oldPoint = newPoint;
		}

		return inside;
	}
	
	public static Vector2? LineSegmentsIntersection(Vector2 p1, Vector2 p2, Vector2 p3, Vector2 p4)
	{
		Vector2? result = null;

		Vector2 a = p2 - p1;
		Vector2 b = p3 - p4;
		Vector2 c = p1 - p3;

		float alphaNumerator = b.y * c.x - b.x * c.y;
		float alphaDenominator = a.y * b.x - a.x * b.y;
		float betaNumerator = a.x * c.y - a.y * c.x;
		float betaDenominator = a.y * b.x - a.x * b.y;

		bool doIntersect = true;

		if (alphaDenominator == 0 || betaDenominator == 0)
		{
			doIntersect = false;
		}
		else
		{

			if (alphaDenominator > 0)
			{
				if (alphaNumerator < 0 || alphaNumerator > alphaDenominator)
				{
					doIntersect = false;
				}
			}
			else if (alphaNumerator > 0 || alphaNumerator < alphaDenominator)
			{
				doIntersect = false;
			}

			if (doIntersect &&  betaDenominator > 0) {
				if (betaNumerator < 0 || betaNumerator > betaDenominator)
				{
					doIntersect = false;
				}
			} else if (betaNumerator > 0 || betaNumerator < betaDenominator)
			{
				doIntersect = false;
			}
		}

		if (doIntersect)
		{

			result = p1 + alphaNumerator * (p2 - p1) / alphaDenominator;
		}
		return result;
	}
}