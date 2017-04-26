using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Drawing;
using System.Drawing.Imaging;

using DispObject;
using DotNetMatrix;

namespace ReinlessLib
{
    public static partial class CReinlessLib
    {
        // Edges From polar croodinate which has most far distance ( Radius = sqrt( imageW^2 + imageH^2);
        public static List<PointF> HC_EDGE_GetPolarEdgeOut(byte[] rawImage, int imageW, int imageH, int nSampling, double fTransX, double fTransY)
        {
            List<PointF> listEdges = new List<PointF>();

            if (rawImage == null) return listEdges;

            int nCenterX = imageW / 2;
            int nCenterY = imageH / 2;

            //int nRadius = Convert.ToInt32(Math.Sqrt((imageW * imageW) + (imageH * imageH)) / 2.0);
            int nRadius = Convert.ToInt32(Math.Sqrt(imageW * imageW + imageH * imageH) /2.0);

            int nAngleStep = 360 / nSampling;
           
            for (int nAngle = 0; nAngle < 360; nAngle += nAngleStep)
            {
                // Inner Contour of deformed circle can not detect when radius is assumed square.
                // +50 additive gap allow to detect deformed circle
                PointF ptEdge = new PointF(-1, -1);
                double fDegree = ((nAngle - 90) * 3.14 / 180);

                for (int moveRadius = nRadius/2; moveRadius < nRadius; moveRadius++)  // START FROM THE CENTER/
                {
                    int x = nCenterX + (int)(moveRadius * Math.Cos(fDegree));
                    int y = nCenterY + (int)(moveRadius * Math.Sin(fDegree));

                    if (x < 0 || y < 0 || x >= imageW || y >= imageH) break;
                    else
                    {
                        if (rawImage[y * imageW + x] == 255)
                        {
                            ptEdge.X = x;
                            ptEdge.Y = y;
                        }
                    }
                }
                if (fTransX != 0 || fTransY != 0)
                {
                    ptEdge = CPoint.OffsetPoint(ptEdge, (float)fTransX, (float)fTransY);
                }
                if (ptEdge.X != -1 && ptEdge.Y != -1)
                {
                    listEdges.Add(ptEdge);
                }
            } // angle loop

            return listEdges;
        } // 161125

        // Edges from normal boundary Edges ( Radius = Max(imageW, imageH )
        public static List<PointF> HC_EDGE_BoundaryEdge(byte[] rawImage, int imageW, int imageH, int nSampling)
        {
            List<PointF> listEdges = new List<PointF>();

            int nCenterX = imageW / 2;
            int nCenterY = imageH / 2;

            //int nRadius = Convert.ToInt32(Math.Sqrt((imageW * imageW) + (imageH * imageH)) / 2.0);
            int nRadius = Math.Max(imageW, imageH);

            byte[] newRaw = new byte[nRadius * 360];

            int nAngleStep = 360 / nSampling;
            for (int nAngle = 0; nAngle < 360; nAngle += nAngleStep)
            {
                int polarY = nRadius;
                // Inner Contour of deformed circle can not detect when radius is assumed square.
                // +50 additive gap allow to detect deformed circle
                PointF ptEdge = new PointF(0, 0);
                for (int moveRadisu = 0; moveRadisu < nRadius; moveRadisu++)
                {
                    double fDegree = 0;
                    fDegree = ((nAngle - 90) * 3.14 / 180);

                    int x = nCenterX + (int)(moveRadisu * Math.Cos(fDegree));
                    int y = nCenterY + (int)(moveRadisu * Math.Sin(fDegree));

                    if (x <= 0 || y <= 0 || x >= imageW || y >= imageH)
                    {
                        continue;
                    }

                    if (rawImage[y * imageW + x] == 255)
                    {
                        ptEdge.X = x;
                        ptEdge.Y = y;
                    }
                } // radius loop
                if (ptEdge.X == 0 || ptEdge.Y == 0) continue; // no matching exception 170105
                listEdges.Add(ptEdge);
            } // angle loop

            return listEdges;
        } // 161125
        public static List<PointF> HC_EDGE_BoundaryEdge(byte[] rawImage, int imageW, int imageH, int nSampling, int nValue) // 170105 for raytopia
        {
            List<PointF> listEdges = new List<PointF>();

            int nCenterX = imageW / 2;
            int nCenterY = imageH / 2;

            //int nRadius = Convert.ToInt32(Math.Sqrt((imageW * imageW) + (imageH * imageH)) / 2.0);
            int nRadius = Math.Max(imageW, imageH);

            byte[] newRaw = new byte[nRadius * 360];

            int nAngleStep = 360 / nSampling;
            for (int nAngle = 0; nAngle < 360; nAngle += nAngleStep)
            {
                int polarY = nRadius;
                // Inner Contour of deformed circle can not detect when radius is assumed square.
                // +50 additive gap allow to detect deformed circle
                PointF ptEdge = new PointF(0, 0);
                for (int moveRadisu = 0; moveRadisu < nRadius; moveRadisu++)
                {
                    double fDegree = 0;
                    fDegree = ((nAngle - 90) * 3.14 / 180);

                    int x = nCenterX + (int)(moveRadisu * Math.Cos(fDegree));
                    int y = nCenterY + (int)(moveRadisu * Math.Sin(fDegree));

                    if (x <= 0 || y <= 0 || x >= imageW || y >= imageH)
                    {
                        continue;
                    }

                    if (rawImage[y * imageW + x] >= nValue)
                    {
                        ptEdge.X = x;
                        ptEdge.Y = y;
                    }
                } // radius loop
                if (ptEdge.X == 0 || ptEdge.Y == 0) continue; // no matching exception 170105
                listEdges.Add(ptEdge);
            } // angle loop

            return listEdges;
        } // 161125

        // variation 1 : rising 0 -> 255
        // variation 0 : falling 255 -> 0
        // dir : to upside = -1
        // dir : to down size = 1
        static public List<PointF> GetHorizontalEdgesBinary(byte[] rawImage, int imageW, int imageH, int variation, int dir = -1)
        {
            List<PointF> listEdges = new List<PointF>();

            int nY = imageH/2;
            if (dir == -1)
            {
                for (int x = 0; x < imageW; x++)
                {
                    for (int y = nY; y > 0; y--)
                    {
                        if (variation == 1 && rawImage[y*imageW+x] == 255)
                        {
                            listEdges.Add(new PointF(x, y));
                            break;
                        }
                        else if( variation == -1 && rawImage[y*imageW+x] == 0 )
                        {
                            listEdges.Add(new PointF(x, y));
                            break;
                        }
                    }
                }
            }
            else if( dir == 1)
            {
                for (int x = 0; x < imageW; x++)
                {
                    for (int y = nY; y < imageH; y++)
                    {
                        if (variation == 1 && rawImage[y * imageW + x] == 255)
                        {
                            listEdges.Add(new PointF(x, y));
                            break;
                        }
                        else if (variation == -1 && rawImage[y * imageW + x] == 0)
                        {
                            listEdges.Add(new PointF(x, y));
                            break;
                        }
                    }
                }
            }
            return listEdges;
        }


       

        /// <summary>
        /// Take the Positional [Left/Top/Right/Bottom] From the arbitrary Point Lists.
        /// 특정 사각 영역안에서 포인트 들이 난재할 경우, 삼각형으로 쪼개서 해당 영역에 포함되는 점들을 구한다.
        /// </summary>
        /// <param name="ptList">arbitrary points </param>
        /// <param name="width"> width of the given region </param>
        /// <param name="height">height of the given region </param>
        /// <returns></returns>
        static public List<PointF> HC_EDGE_GetPointsTOP(List<PointF> ptList, int width, int height)
        {
            List<PointF> ptListTop = new List<PointF>();

            PointF ptA = new PointF(0, 0);
            PointF ptB = new PointF(width, 0);
            PointF ptC = new PointF(width / 2, height / 2);

            for (int i = 0; i < ptList.Count; i++)
            {
                if (HC_EDGE_IsIntersectTriangleAndPoint(ptList.ElementAt(i), ptA, ptB, ptC) == true)
                {
                    ptListTop.Add(ptList.ElementAt(i));
                }
            }
            return ptListTop;
        }
        static public List<PointF> HC_EDGE_GetPointsBTM(List<PointF> ptList, int width, int height)
        {
            List<PointF> ptListBTM = new List<PointF>();

            PointF ptA = new PointF(width / 2, height / 2);
            PointF ptB = new PointF(width, height);
            PointF ptC = new PointF(0, height);

            for (int i = 0; i < ptList.Count; i++)
            {
                if (HC_EDGE_IsIntersectTriangleAndPoint(ptList.ElementAt(i), ptA, ptB, ptC) == true)
                {
                    ptListBTM.Add(ptList.ElementAt(i));
                }
            }
            return ptListBTM;
        }
        static public List<PointF> HC_EDGE_GetPointsLFT(List<PointF> ptList, int width, int height)
        {
            List<PointF> ptListLFT = new List<PointF>();

            PointF ptA = new PointF(0, 0);
            PointF ptB = new PointF(width / 2, height / 2);
            PointF ptC = new PointF(0, height);

            for (int i = 0; i < ptList.Count; i++)
            {
                if (HC_EDGE_IsIntersectTriangleAndPoint(ptList.ElementAt(i), ptA, ptB, ptC) == true)
                {
                    ptListLFT.Add(ptList.ElementAt(i));
                }
            }
            return ptListLFT;
        }
        static public List<PointF> HC_EDGE_GetPointsRHT(List<PointF> ptList, int width, int height)
        {
            List<PointF> ptListRHT = new List<PointF>();

            PointF ptA = new PointF(width, 0);
            PointF ptB = new PointF(width, height);
            PointF ptC = new PointF(width / 2, height / 2);

            for (int i = 0; i < ptList.Count; i++)
            {
                if (HC_EDGE_IsIntersectTriangleAndPoint(ptList.ElementAt(i), ptA, ptB, ptC) == true)
                {
                    ptListRHT.Add(ptList.ElementAt(i));
                }
            }
            return ptListRHT;
        }
        static public bool HC_EDGE_IsIntersectTriangleAndPoint(PointF ptS, PointF ptA, PointF ptB, PointF ptC)
        {
            double diffSA_X = ptS.X - ptA.X;
            double diffSA_Y = ptS.Y - ptA.Y;

            bool s_ab = (ptB.X - ptA.X) * diffSA_Y - (ptB.Y - ptA.Y) * diffSA_X > 0;

            if ((ptC.X - ptA.X) * diffSA_Y - (ptC.Y - ptA.Y) * diffSA_X > 0 == s_ab) return false;

            if ((ptC.X - ptB.X) * (ptS.Y - ptB.Y) - (ptC.Y - ptB.Y) * (ptS.X - ptB.X) > 0 != s_ab) return false;

            return true;
        }

        // 최장 거리 라인 추출 -->  최장 거리 라인 대비 대각선 라인 추출
        // 사각형일경우 최장축 기준 대각이 거의 유사해야 되니까.
        static public void /****/HC_EDGE_GetCrossLine(List<PointF> ptList, ref CLine l1, ref CLine l2)
        {
            CLine line1 = new CLine();
            CLine line2 = new CLine();

            Parallel.Invoke(()=>
            {
                int nHalf = ptList.Count / 2;

                List<CLine> listLines = new List<CLine>();

                var listLinePair = new List<KeyValuePair<double, int>>();

                for (int nIndex = 0; nIndex < nHalf; nIndex++)
                {
                    CLine crossLine = new CLine(ptList.ElementAt(nIndex), ptList.ElementAt(nIndex + nHalf));
                    listLines.Add(crossLine);
                    listLinePair.Add(new KeyValuePair<double, int>(crossLine.LENGTH, nIndex));
                }

                if (listLinePair.Count < 2) return;

                listLinePair.Sort(SortDoubleInt);

                // get the most longest line 
                int index1 = listLinePair.ElementAt(0).Value;
                line1 = listLines.ElementAt(index1);

                // get the digonal angle with sign compensation 
                int nPos = 0;
                if (index1 < 90)
                {
                    nPos = 90 + index1;
                }
                else if (index1 >= 90)
                {
                    nPos = index1 - 90;
                }

                // search the longest digonal line  within +/- 10 degree
                int nStart = nPos - 10; nStart = nStart < 0 ? 0 : nStart;
                int nEnd = nPos + 10; nEnd = nEnd > 180 ? 180 : nEnd;

                double max = 0;
                int maxPos = 0;
                for (int loop = nStart; loop < nEnd; loop++)
                {
                    line2 = listLines.ElementAt(loop);

                    if (line2.LENGTH > max)
                    {
                        max = line2.LENGTH;
                        maxPos = loop;
                    }
                }

                line2 = listLines.ElementAt(maxPos);
            //}
            });
            l1 = line1;
            l2 = line2;
        }

        static int SortDoubleInt(KeyValuePair<double, int> a, KeyValuePair<double, int> b)
        {
            return b.Key.CompareTo(a.Key);
        }

        static public CLine HC_EDGE_FitLineHor(List<PointF> listPoints, int w, int h)
        {
            double meanX = 0;
            double meanY = 0;

            for (int i = 0; i < listPoints.Count; i++)
            {
                meanX += listPoints.ElementAt(i).X;
                meanY += listPoints.ElementAt(i).Y;
            }
            meanX /= listPoints.Count;
            meanY /= listPoints.Count;

            double numer = 0;
	        double denom = 0;

            for( int i = 0; i < listPoints.Count; i++)
            {
                numer += (listPoints.ElementAt(i).X - meanX ) * (listPoints.ElementAt(i).Y - meanY);
                denom += (listPoints.ElementAt(i).X-meanX)* (listPoints.ElementAt(i).X-meanX);
            }
            double m = numer / denom;
	        double b = meanY - (m*meanX);

            CLine line = new CLine();
            line.P1 = (new PointF(0, (float)((m * line.P1.X) + b)));
            line.P2 = (new PointF(w, (float)((m * line.P2.X) + b)));

            return line;
	    }
        static public CLine HC_EDGE_FitLineVer(List<PointF> listPoints, int w, int h)
        {
            List<PointF> ptTemp = new List<PointF>();

            ptTemp = HC_EDGE_GetRotatedPoints(listPoints, 90, w, h);

            double meanX = 0;
            double meanY = 0;

            for (int i = 0; i < ptTemp.Count; i++)
            {
                meanX += ptTemp.ElementAt(i).X;
                meanY += ptTemp.ElementAt(i).Y;
            }
            meanX /= ptTemp.Count;
            meanY /= ptTemp.Count;

            double numer = 0;
            double denom = 0;

            for (int i = 0; i < ptTemp.Count; i++)
            {
                numer += (ptTemp.ElementAt(i).X - meanX) * (ptTemp.ElementAt(i).Y - meanY);
                denom += (ptTemp.ElementAt(i).X - meanX) * (ptTemp.ElementAt(i).X - meanX);
            }
            double m = numer / denom;
            double b = meanY - (m * meanX);

            CLine line = new CLine();

            line.P1 = new PointF(0, (float)((line.P1.Y - b) / m));
            line.P2 = new PointF((float)((line.P2.Y - b) / m), h);

            return line.RotateLinebyPoints(-90, w/2, h/2);

        }

        /// <summary>
        /// rotate points by input angle based on center x,y
        /// </summary>
        /// <param name="listPoints"></param>
        /// <param name="fAngle"></param>
        /// <param name="w">to calculate center x</param>
        /// <param name="h">to calculate center y</param>
        /// <returns></returns>
        static public List<PointF> HC_EDGE_GetRotatedPoints(List<PointF> listPoints, double fAngle, int w, int h)
        {
            List<PointF> ptRotated = new List<PointF>();

            double fDegree = fAngle * Math.PI / 180.0;

            for (int i = 0; i < listPoints.Count; i++)
            {
                Point pt = _GetRotatePos(listPoints.ElementAt(i).X, listPoints.ElementAt(i).Y, fDegree, w / 2, h / 2);

                ptRotated.Add(new PointF(pt.X, pt.Y));
            }
            return ptRotated;
        }

        public static void HC_FIT_Circle(List<PointF> list, ref PointF ptCenter, ref double radius)
        {
            double sx = 0.0, sy = 0.0;
            double sx2 = 0.0, sy2 = 0.0, sxy = 0.0;
            double sx3 = 0.0, sy3 = 0.0, sx2y = 0.0, sxy2 = 0.0;

            /* compute summations */
            for (int k = 0; k < list.Count; k++)
            {
                double x = list.ElementAt(k).X;
                double y = list.ElementAt(k).Y;

                double xx = x * x;
                double yy = y * y;

                sx = sx + x;
                sy = sy + y;
                sx2 = sx2 + xx;
                sy2 = sy2 + yy;
                sxy = sxy + x * y;
                sx3 = sx3 + x * xx;
                sy3 = sy3 + y * yy;
                sx2y = sx2y + xx * y;
                sxy2 = sxy2 + yy * x;
            }
            /* compute a's,b's,c's */
            double a1 = 2.0 * (sx * sx - sx2 * list.Count);
            double a2 = 2.0 * (sx * sy - sxy * list.Count);
            double b1 = a2;
            double b2 = 2.0 * (sy * sy - sy2 * list.Count);
            double c1 = sx2 * sx - sx3 * list.Count + sx * sy2 - sxy2 * list.Count;
            double c2 = sx2 * sy - sy3 * list.Count + sy * sy2 - sx2y * list.Count;

            double det = a1 * b2 - a2 * b1;
            if (Math.Abs(det) < 0.0001)
            {                /*collinear한 경우임;*/
                return;
            }

            /* floating value  center */
            double cx = (c1 * b2 - c2 * b1) / det;
            double cy = (a1 * c2 - a2 * c1) / det;

            /* compute radius squared */
            double radsq = (sx2 - 2 * sx * cx + cx * cx * list.Count + sy2 - 2 * sy * cy + cy * cy * list.Count) / list.Count;
            radius = Math.Sqrt(radsq);
            /* integer value center */
            ptCenter.X = Convert.ToSingle(cx + 0.5);
            ptCenter.Y = Convert.ToSingle(cy + 0.5);

            return;
        }
        public static List<PointF> HC_FIT_Ellipse(List<PointF> ptListTarget, int density, ref PointF ptCenter)
        {
            double A = 0, B = 0;
            double cos_phi = 0, sin_phi = 0;

            double ptCX = 0;
            double ptCY = 0;
            HC_FIT_EllipseParamSet(ptListTarget, out ptCX, out ptCY, out A, out B, out cos_phi, out sin_phi);

            List<PointF> ptContourList = HC_FIT_EllipseGenContour(density, ptCX, ptCY, A, B, cos_phi, sin_phi);

            double fAVG_X = 0;
            double fAVG_Y = 0;

            foreach (PointF pt in ptContourList)
            {
                fAVG_X += pt.X; fAVG_Y += pt.Y;
            }
            fAVG_X /= ptContourList.Count;
            fAVG_Y /= ptContourList.Count;

            ptCenter = new PointF((float)fAVG_X, (float)fAVG_Y);

            return ptContourList;
        }
        public static List<PointF> HC_FIT_EllipseGenContour(int nSamplingDensity, double CX, double CY, double A, double B, double cos_phi, double sin_phi)
        {
            double fPitch = (2 * Math.PI) / nSamplingDensity;

            double[] contArrX = new double[nSamplingDensity];
            double[] contArrY = new double[nSamplingDensity];

            Parallel.For(0, nSamplingDensity, i =>
            {
                contArrX[i] = CX + (A * Math.Cos(fPitch * i));
                contArrY[i] = CY + (B * Math.Sin(fPitch * i));
            });

            double[][] arrContour = { contArrX, contArrY };
            DotNetMatrix.Matrix mtrpreContour = new DotNetMatrix.Matrix(arrContour);

            double[][] arrSC = { new double[] { cos_phi, sin_phi }, new double[] { -sin_phi, cos_phi } };
            DotNetMatrix.Matrix mtrSC = new DotNetMatrix.Matrix(arrSC);

            DotNetMatrix.Matrix mtrContour = mtrSC.Multiply(mtrpreContour);

            float x = 0;
            float y = 0;

            List<PointF> ptList = new List<PointF>();

            for (int i = 0; i < nSamplingDensity; i++)
            {
                x = (float)mtrContour.GetElement(0, i);
                y = (float)mtrContour.GetElement(1, i);

                ptList.Add(new PointF(x, y));
            }

            return ptList;
        }
        public static void HC_FIT_EllipseParamSet(List<PointF> ptListTarget, out double ptCX, out double ptCY, out double A, out double B, out double cos_phi, out double sin_phi)
        {
            if (ptListTarget.Count == 0)
            {
                A = B = 0;
                cos_phi = sin_phi = 0;
                ptCX = ptCY = 0;
                return;
            }
            int nDataCount = ptListTarget.Count;

            float[] CroodList_X = ptListTarget.Select(element => element.X).ToArray();
            float[] CroodList_Y = ptListTarget.Select(element => element.Y).ToArray();

            // required output : a, b, sign, cos, center x, center y
            A = B = 0; cos_phi = sin_phi = 0; ptCX = ptCY = 0;

            double meanX = CroodList_X.Average();
            double meanY = CroodList_Y.Average();

            double[] x1 = new double[nDataCount];
            double[] y1 = new double[nDataCount];
            double[] xx = new double[nDataCount];
            double[] yy = new double[nDataCount];
            double[] xy = new double[nDataCount];

            for (int i = 0; i < nDataCount; i++)
            {
                x1[i] = CroodList_X[i] - meanX;
                y1[i] = CroodList_Y[i] - meanY;
            }

            for (int i = 0; i < nDataCount; i++)
            {
                xx[i] = x1[i] * x1[i];
                yy[i] = y1[i] * y1[i];
                xy[i] = x1[i] * y1[i];
            }

            double[][] arrFittingData = { xx, xy, yy, x1, y1 };
            double[] arrFittingDataSum = new double[5];

            arrFittingDataSum[0] = xx.Sum();
            arrFittingDataSum[1] = xy.Sum();
            arrFittingDataSum[2] = yy.Sum();
            arrFittingDataSum[3] = x1.Sum();
            arrFittingDataSum[4] = y1.Sum();

            // make nemerator *********************************************************************
            // sum ( fitting data 1X5 matrix )

            DotNetMatrix.Matrix mtrNumerator = new DotNetMatrix.Matrix(arrFittingDataSum, 1);

            // make denominator  ******************************************************************
            // fitting data' * fitting data   by transpose & multiply

            DotNetMatrix.Matrix mtrFittingData = new DotNetMatrix.Matrix(arrFittingData);
            DotNetMatrix.Matrix mtrTrans = mtrFittingData.Transpose();
            DotNetMatrix.Matrix mtrDenominator = mtrFittingData.Multiply(mtrTrans);
            DotNetMatrix.Matrix mtrInvertDenom = mtrDenominator.Inverse();

            // calcuate parameters ****************************************************************
            DotNetMatrix.Matrix res = mtrNumerator.Multiply(mtrInvertDenom);

            A = res.GetElement(0, 0);
            B = res.GetElement(0, 1);

            double C = res.GetElement(0, 2);
            double D = res.GetElement(0, 3);
            double E = res.GetElement(0, 4);

            double oriental_Rad = 0;

            // fitting check value 
            double checkValue = Math.Min(Math.Abs(B / A), Math.Abs(B / C));

            if (checkValue > 0.0001) // pass case 
            {
                oriental_Rad = 1.0 / 2.0 * Math.Atan(B / (C - A));

                cos_phi = Math.Cos(oriental_Rad);
                sin_phi = Math.Sin(oriental_Rad);

                A = (A * cos_phi * cos_phi) - (B * cos_phi * sin_phi) + (C * sin_phi * sin_phi);
                B = 0;
                C = (A * sin_phi * sin_phi) + (B * cos_phi * sin_phi) + (C * cos_phi * cos_phi);
                D = (D * cos_phi - E * sin_phi);
                E = (D * sin_phi + E * cos_phi);

                meanX = cos_phi * meanX - sin_phi * meanY;
                meanY = sin_phi * meanX + cos_phi * meanY;
            }
            else // false case 
            {
                oriental_Rad = 0;
                cos_phi = Math.Cos(oriental_Rad);
                sin_phi = Math.Sin(oriental_Rad);
            }

            double Status = A * C;
            double F = 0;

            /***/
            if (Status == 0) Console.Write("Parabola Found\n");
            else if (Status < 0) Console.Write("Hyperbola Found\n");
            else if (Status > 0)
            {
                if (A < 0)
                {
                    A *= -1; C *= -1; D *= -1; E *= -1;
                }
                else
                {
                    ptCX = (float)(meanX - D / 2 / A);
                    ptCY = (float)(meanY - E / 2 / C);

                    F = 1 + (D * D) / (4 * A) + (E * E) / (4 * C);

                    A = Math.Sqrt(F / A);
                    B = Math.Sqrt(F / C);

                    #region MyRegion
                    // no meaning full variables for axis
                    //double long_Axis = 2 * Math.Max(a, b);
                    //double shor_Axis = 2 * Math.Min(a, b);

                    // Circle center calculation
                    //double[] arrXY = new double[] { CX, CY };

                    //Matrix mtrPosIN = new Matrix(arrXY, 1);
                    //Matrix mtrPosRes = mtrPosIN.Multiply(mtrSC);

                    //X0_in = (X0 * cos_phi )+  (Y0 * sin_phi);
                    //Y0_in = (X0 * -sin_phi) + (Y0 * cos_phi);

                    //ptCenter.X = Convert.ToInt32(mtrPosRes.GetElement(0, 0));
                    //ptCenter.Y = Convert.ToInt32(mtrPosRes.GetElement(0, 1));
                    #endregion

                    #region cross line calculation - no need

                    //double[][] verLine = { new double[] { X0, X0 }, new double[] { Y0 + (-B), Y0 + B } };
                    //double[][] horLine = { new double[] { X0 + (-A), X0 + A }, new double[] { Y0, Y0 } };

                    //Matrix mtrVerLine = new Matrix(verLine);
                    //Matrix mtrHorLine = new Matrix(horLine);

                    //Matrix mtrFinalLineVer = mtrSC.Multiply(mtrVerLine);
                    //Matrix mtrFinalLineHor = mtrSC.Multiply(mtrHorLine);

                    //Point pt_v1 = new Point(); Point pt_v2 = new Point();
                    //Point pt_h1 = new Point(); Point pt_h2 = new Point();

                    //pt_v1.X = (int)(mtrFinalLineVer.GetElement(0, 0));
                    //pt_v1.Y = (int)(mtrFinalLineVer.GetElement(1, 0));
                    //pt_v2.X = (int)(mtrFinalLineVer.GetElement(0, 1));
                    //pt_v2.Y = (int)(mtrFinalLineVer.GetElement(1, 1));

                    //pt_h1.X = (int)(mtrFinalLineHor.GetElement(0, 0));
                    //pt_h1.Y = (int)(mtrFinalLineHor.GetElement(1, 0));
                    //pt_h2.X = (int)(mtrFinalLineHor.GetElement(0, 1));
                    //pt_h2.Y = (int)(mtrFinalLineHor.GetElement(1, 1));

                    #endregion

                }
            }
        }

        /// <summary>
        /// 특정 형태 영역(Elliptical)에 포함된 영역의 점들을 털어내서 구한다. 
        /// </summary>
        /// <param name="rc"></param>
        /// <param name="list"></param>
        /// <returns></returns>
        public static List<PointF> GetFilteredEllipsePoints(RectangleF rc, List<PointF> list)
        {
            System.Drawing.Drawing2D.GraphicsPath myPath = new System.Drawing.Drawing2D.GraphicsPath();
            myPath.AddEllipse(rc);

            List<PointF> listTemp = new List<PointF>();

            for (int i = 0; i < list.Count; i++)
            {
                PointF pt = list.ElementAt(i);
                if (myPath.IsVisible(pt) == true) { listTemp.Add(pt); }
            }
            return listTemp;
        }

        public static double Get2ndDerivativeLine_MaxPos(byte[] line)
        {
            double[] buff_1st = new double[line.Length - 1];
            double[] buff_2nd = new double[line.Length - 2];

            for (int nIndex = 0; nIndex < line.Length - 1; nIndex++) { buff_1st[nIndex] = line[nIndex + 1] - line[nIndex]; }
            for (int nIndex = 0; nIndex < line.Length - 2; nIndex++) { buff_2nd[nIndex] = buff_1st[nIndex + 1] - buff_1st[nIndex]; }

            double fMax = buff_2nd.Max();
            int/**/nPos = Array.IndexOf(buff_2nd, fMax);

            double fSubPixel = 0;
            try { fSubPixel = _GetSubPixel(buff_2nd.ElementAt(nPos - 1), buff_2nd.ElementAt(nPos), buff_2nd.ElementAt(nPos + 1)); }
            catch { }

            return fSubPixel + nPos;
        }
        public static double Get2ndDerivativeLine_MinPos(byte[] line)
        {
            double[] buff_1st = new double[line.Length - 1];
            double[] buff_2nd = new double[line.Length - 2];

            for (int nIndex = 0; nIndex < line.Length - 1; nIndex++) { buff_1st[nIndex] = line[nIndex + 1] - line[nIndex]; }
            for (int nIndex = 0; nIndex < line.Length - 2; nIndex++) { buff_2nd[nIndex] = buff_1st[nIndex + 1] - buff_1st[nIndex]; }

            double fMin = buff_2nd.Min();
            int/**/nPos = Array.IndexOf(buff_2nd, fMin);

            double fSubPixel = 0;
            try { fSubPixel = _GetSubPixel(buff_2nd.ElementAt(nPos - 1), buff_2nd.ElementAt(nPos), buff_2nd.ElementAt(nPos + 1)); }
            catch { }

            return fSubPixel + nPos;
        }

        public static List<PointF> Get2ndDerivativeList_HorMax(byte[] rawImage, int imageW, int imageH)
        {
            double[] buff_1st = new double[imageW - 1];
            double[] buff_2nd = new double[imageW - 2];
            double fSubPixel = 0;

            List<PointF> list = new List<PointF>();

            for (int y = 0; y < imageH; y++)
            {
                Array.Clear(buff_1st, 0, buff_1st.Length);
                Array.Clear(buff_2nd, 0, buff_2nd.Length);
                fSubPixel = 0;

                for (int x = 0; x < imageW - 1; x++)
                {
                    buff_1st[x] = rawImage[y * imageW + x + 1] - rawImage[y * imageW + x];
                }

                for (int x = 0; x < imageW - 2; x++)
                {
                    buff_2nd[x] = buff_1st[x + 1] - buff_1st[x];
                }
                double fMax = buff_2nd.Max();
                int nPosMax = Array.IndexOf(buff_2nd, fMax);
                try { fSubPixel = _GetSubPixel(buff_2nd.ElementAt(nPosMax - 1), buff_2nd.ElementAt(nPosMax), buff_2nd.ElementAt(nPosMax + 1)); }
                catch { }

                list.Add(new PointF((float)(nPosMax + fSubPixel), y));

            }
            return list;
        }
        public static List<PointF> Get2ndDerivativeList_HorMin(byte[] rawImage, int imageW, int imageH)
        {
            double[] buff_1st = new double[imageW - 1];
            double[] buff_2nd = new double[imageW - 2];
            double fSubPixel = 0;

            List<PointF> list = new List<PointF>();

            for (int y = 0; y < imageH; y++)
            {
                Array.Clear(buff_1st, 0, buff_1st.Length);
                Array.Clear(buff_2nd, 0, buff_2nd.Length);
                fSubPixel = 0;

                for (int x = 0; x < imageW - 1; x++)
                {
                    buff_1st[x] = rawImage[y * imageW + x + 1] - rawImage[y * imageW + x];
                }

                for (int x = 0; x < imageW - 2; x++)
                {
                    buff_2nd[x] = buff_1st[x + 1] - buff_1st[x];
                }
                double fMin = buff_2nd.Min();
                int nPosMin = Array.IndexOf(buff_2nd, fMin);
                try { fSubPixel = _GetSubPixel(buff_2nd.ElementAt(nPosMin - 1), buff_2nd.ElementAt(nPosMin), buff_2nd.ElementAt(nPosMin + 1)); }
                catch { }

                list.Add(new PointF((float)(nPosMin + fSubPixel), y));

            }
            return list;
        }
        public static double _GetSubPixel(double pa, double pb, double pc)
        {
            // simple quadratic interpolation
            return 0.5 * (pa - pc) / (pa - (2 * pb) + pc);
        }

        // get the X difference from the list points    161129
        static public double HC_EDGE_GetDifferenceX(List<PointF> listPoints)
        {
            double[] arrPosY = new double[listPoints.Count];

            for (int i = 0; i < listPoints.Count; i++)
            {
                arrPosY[i] = listPoints.ElementAt(i).X;
            }

            return arrPosY.Max() - arrPosY.Min();

        }
        // get the Y difference from the list points 161129
        static public double HC_EDGE_GetDifferenceY(List<PointF> listPoints)
        {
            PointF[] arrPoints = listPoints.ToArray();
            double[] arrPosY = new double[arrPoints.Length];

            for (int i = 0; i < arrPoints.Length; i++)
            {
                arrPosY[i] = arrPoints[i].Y;
            }

            return arrPosY.Max() - arrPosY.Min();
        }


        public static double GetZeroCrossingPt(byte [] buff, int sign, double mag )
        {
            int [] differFST = new int[buff.Length];
            int [] differSCD = new int[buff.Length];

            for (int i = 0; i < buff.Length - 1; i++){differFST[i] = buff[i + 1] - buff[i];}
            for (int i = 0; i < buff.Length - 2; i++){differSCD[i] = differFST[i + 1] - differFST[i];}

             double ptEdge = 0.0;

             for (int j = 0; j < buff.Length - 2; j++)
             {
                 if (((sign == -1) && (differFST[j + 1] < sign * mag)) /*|| 
                     ((sign == +1) && (differFST[j + 1] > sign * mag)) */)
                 {
                     if ((differSCD[j] * differSCD[j + 1] < 0))
                     {
                         double subpixel_range = (double)(-differSCD[j]) / (differSCD[j + 1] - differSCD[j]);
                         ptEdge = (j + 1) + subpixel_range;
                         break;
                     }
                     else if (j > 0 && differSCD[j] == 0)
                     {
                         if (differSCD[j - 1] * differSCD[j + 1] < 0)
                         {
                             ptEdge = j + 1;
                             break;
                         }
                     }
                 }
             }

             if(Math.Floor(ptEdge) == Math.Floor(0.0))
             {
                 for (int j = 0; j < buff.Length - 2; j++)
                 {
                     if (((sign == -1) && (differFST[j + 1] > -sign * mag)) /*|| 
                         ((sign == +1) && (differFST[j + 1] < -sign * mag))*/)
                     {
                         if ((differSCD[j] * differSCD[j + 1] < 0))
                         {
                             double subpixel_range = (double)(-differSCD[j]) / (differSCD[j + 1] - differSCD[j]);
                             ptEdge = (j + 1) + subpixel_range;
                             break;
                         }
                         else if (j > 0 && differSCD[j] == 0)
                         {
                             if (differSCD[j - 1] * differSCD[j + 1] < 0)
                             {
                                 ptEdge = j + 1;
                                 break;
                             }
                         }
                     }
                 }
             }
             return ptEdge;
         }


        public static float[] GetElements_X(List<PointF> list)
        {
            float[] arrX = (float[])list.Select(element => element.X).ToArray();
            return arrX;
        }
        public static float[] GetElements_Y(List<PointF> list)
        {
            float[] arrY = (float[])list.Select(element => element.Y).ToArray();
            return arrY;
        }


        // fitering functions for the horizontal and vertical rectangle boundary 170412 
        public static List<PointF> GetList_FilteredBy_TY_BY(List<PointF> list, double fTY, double fBY, double fThreshold)
        {
            List<PointF> listBuff = new List<PointF>();

            for (int i = 0; i < list.Count; i++)
            {
                PointF pt = list.ElementAt(i);
                if (Math.Abs(pt.Y - fTY) < fThreshold || Math.Abs(pt.Y - fBY) < fThreshold) continue;
                listBuff.Add(pt);
            }
            return listBuff;
        }
        public static List<PointF> GetList_FilteredBy_LX_RX(List<PointF> list, double fLX, double fRX, double fThreshold)
        {
            List<PointF> listBuff = new List<PointF>();

            for (int i = 0; i < list.Count; i++)
            {
                PointF pt = list.ElementAt(i);
                if (Math.Abs(pt.X - fLX) < fThreshold || Math.Abs(pt.X - fRX) < fThreshold) continue;
                listBuff.Add(pt);
            }
            return listBuff;
        }
        public static List<PointF> GetList_FilterBy_MajorDistance(List<PointF> list, bool bAxisX, double fDistance)
        {
            List<PointF> listBuff = new List<PointF>();

            PointF[] arrPoints = list.ToArray();

            int nMax = 0;
            int nMaxPos = 0;

            if (bAxisX == true)
            {

                int[] arrX = arrPoints.Select(element => (int)element.X).ToArray();
                nMax = arrX.Max() + 1;
                int[] nHisto = new int[nMax];

                for (int i = 0; i < arrX.Length; i++)
                {
                    nHisto[arrX[i]]++;
                }
                nMax = nHisto.Max();
                nMaxPos = Array.IndexOf(nHisto, nMax);

                for (int i = 0; i < list.Count; i++)
                {
                    PointF pt = list.ElementAt(i);
                    if (Math.Abs(pt.X - nMaxPos) < fDistance)
                    {
                        listBuff.Add(pt);
                    }
                    else
                    {
                        listBuff.Add(new PointF(nMaxPos, pt.Y));
                    }
                }
            }
            else if (bAxisX == false)
            {
                int[] arrY = arrPoints.Select(element => (int)element.Y).ToArray();
                nMax = arrY.Max() + 1;
                int[] nHisto = new int[nMax];

                for (int i = 0; i < arrY.Length; i++)
                {
                    nHisto[arrY[i]]++;
                }
                nMax = nHisto.Max();
                nMaxPos = Array.IndexOf(nHisto, nMax);

                for (int i = 0; i < list.Count; i++)
                {
                    PointF pt = list.ElementAt(i);
                    if (Math.Abs(pt.Y - nMaxPos) < fDistance)
                    {
                        listBuff.Add(pt);
                    }
                    else
                    {
                        listBuff.Add(new PointF(pt.X, nMaxPos));
                    }
                }
            }
            return listBuff;
        }

        //********************************************************************************************************************************
        // POINT ENFORCING FUNCTIONS 
        //********************************************************************************************************************************
        // croodinate points fixation function for the axis  170412 
        public static List<PointF> ReplacePointList_Absolute_X(RectangleF rc, float x)
        {
            List<PointF> list = new List<PointF>();

            int nHead = (int)rc.Y;
            int nTail = (int)rc.Y + (int)rc.Height;

            for (int y = nHead; y < nTail; y++)
            {
                list.Add(new PointF(x, y));
            }
            return list;
        }
        public static List<PointF> ReplacePointList_Absolute_Y(RectangleF rc, float y)
        {
            List<PointF> list = new List<PointF>();

            int nHead = (int)rc.X;
            int nTail = (int)rc.X + (int)rc.Width;

            for (int x = nHead; x < nTail; x++)
            {
                list.Add(new PointF(x, y));
            }
            return list;
        }

        //CPoint CvvImage::GetIntersecPoint(CPoint H1, CPoint H2, CPoint V1, CPoint V2)
        //{
        //    double fDx1 = H2.x - H1.x;
        //    double fDy1 = H2.y - H1.y;

        //    double fDx2 = V2.x - V1.x;
        //    double fDy2 = V2.y - V1.y;

        //    double b_dot_d_perp = fDx1 * fDy2 - fDy1 * fDx2;


        //    CPoint ptIntersect = CPoint(0, 0);

        //    if (b_dot_d_perp == 0)
        //    {
        //        return ptIntersect;
        //    }

        //    double cx = V1.x - H1.x;
        //    double cy = V1.y - H1.y;

        //    double t = (cx * fDy2 - cy * fDx2) / b_dot_d_perp;

        //    ptIntersect.x = (int)(H1.x + t * fDx1);
        //    ptIntersect.y = (int)(H1.y + t * fDy1);

        //    return ptIntersect;
        //}

        /** 
 
         * */


    }
}
