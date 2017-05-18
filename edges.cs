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
        static public List<PointF> HC_EDGE_GetHorizontalEdgesBinary(byte[] rawImage, int imageW, int imageH, int variation, int dir = -1)
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


        //******************************************************************************************
        // Point selection
        //******************************************************************************************

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
                if (_IsIntersectTriangleAndPoint(ptList.ElementAt(i), ptA, ptB, ptC) == true)
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
                if (_IsIntersectTriangleAndPoint(ptList.ElementAt(i), ptA, ptB, ptC) == true)
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
                if (_IsIntersectTriangleAndPoint(ptList.ElementAt(i), ptA, ptB, ptC) == true)
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
                if (_IsIntersectTriangleAndPoint(ptList.ElementAt(i), ptA, ptB, ptC) == true)
                {
                    ptListRHT.Add(ptList.ElementAt(i));
                }
            }
            return ptListRHT;
        }
        static public bool _IsIntersectTriangleAndPoint(PointF ptS, PointF ptA, PointF ptB, PointF ptC)
        {
            double diffSA_X = ptS.X - ptA.X;
            double diffSA_Y = ptS.Y - ptA.Y;

            bool s_ab = (ptB.X - ptA.X) * diffSA_Y - (ptB.Y - ptA.Y) * diffSA_X > 0;

            if ((ptC.X - ptA.X) * diffSA_Y - (ptC.Y - ptA.Y) * diffSA_X > 0 == s_ab) return false;

            if ((ptC.X - ptB.X) * (ptS.Y - ptB.Y) - (ptC.Y - ptB.Y) * (ptS.X - ptB.X) > 0 != s_ab) return false;

            return true;
        }

        /// <summary>
        /// 특정 형태 영역(Elliptical)에 포함된 영역의 점들을 털어내서 구한다. 
        /// </summary>
        /// <param name="rc"></param>
        /// <param name="list"></param>
        /// <returns></returns>
        public static List<PointF> HC_EDGE_GetFilteredEllipsePoints(RectangleF rc, List<PointF> list)
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


        //******************************************************************************************
        // Fitting Functions 
        //******************************************************************************************

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

        public static void /******/HC_FIT_Circle(List<PointF> list, ref PointF ptCenter, ref double radius)
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
            _EllipseParamSet(ptListTarget, out ptCX, out ptCY, out A, out B, out cos_phi, out sin_phi);

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
        public static void _EllipseParamSet(List<PointF> ptListTarget, out double ptCX, out double ptCY, out double A, out double B, out double cos_phi, out double sin_phi)
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

        //*******************************************************************************************
        // 2nd Derivative based Edge detection functions 
        //*******************************************************************************************
        
        public static double[] HC_EDGE_Get2ndDerivativeArrayFromLineBuff(double[] fLineBuff)
        {
            double[] arr1st = new double[fLineBuff.Length - 1];
            double[] arr2nd = new double[fLineBuff.Length - 2];

            for (int i = 0; i < fLineBuff.Length - 1; i++)
            {
                arr1st[i] = fLineBuff[i + 1] - fLineBuff[i];
            }
            for (int i = 0; i < arr1st.Length - 1 + 0; i++)
            {
                arr2nd[i] = arr1st[i + 1] - arr1st[i];
            }
            return arr2nd;
        }
        public static double HC_EDGE_Get2ndDerivativeLine_MaxPos(byte[] line)
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
        public static double HC_EDGE_Get2ndDerivativeLine_MinPos(byte[] line)
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

        // 대각선같이 기준선이 있을때 측정 축에 대하여 만족하는 2차미분 포지션 추출
        // 후처리로, 특정 기준선(사선) 개별 포인트에 해당 위치를 더해서, 사선엣지를 구할때 쓴다. 
        public static List<PointF> HC_EDGE_Get2ndDerivativeList_HorMax(byte[] rawImage, int imageW, int imageH)
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
        public static List<PointF> HC_EDGE_Get2ndDerivativeList_HorMin(byte[] rawImage, int imageW, int imageH)
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

        // 3개의 값을 기준으로 (전,중,후) 가중치 값을 계산, 2차 미분값 위치가 정해진 경우, 미세보정 용도
        public static double _GetSubPixel(double pa, double pb, double pc)
        {
            // simple quadratic interpolation
            return 0.5 * (pa - pc) / (pa - (2 * pb) + pc);
        }

        // relative function 1: GetPojection[POS]_For[DIR]_Derivative { Get Projection for each direction }
        // relative function 2: ReplaceDerivative_[DIR]_Average  { Replace by representative Position and make point list} 170419 
        // 2차 미분을 위해서, 영역과 방향에 따라, 평균 projection 값 생성해서 축에 대한 단일 2차 미분 array 반환 
        public static double[]/***/HC_EDGE_GetProjection_TOP_For_Hor_Derivative(byte[] rawImage, int imageW, int imageH, RectangleF rc, int nDir)
        {
            const int DIR_INFALL = 0;      // Direction = To Inside Falling
            const int DIR_INRISE = 1;      // Direction = to Inside Rising
            const int DIR_EXFALL = 2;      // Direction = to outside Falling
            const int DIR_EXRISE = 3;      // direction = to Outside Rising

            int sx = (int)rc.X;
            int ex = (int)rc.Width + sx;
            int sy = (int)rc.Y;
            int ey = (int)rc.Height + sy;

            double[] buff_Org = new double[(int)rc.Height + 2];

            if (nDir == DIR_INFALL || nDir == DIR_INRISE)
            {
                for (int y = sy; y < ey + 2; y++)
                {
                    for (int x = sx; x < ex; x++)
                    {
                        buff_Org[y - sy] += rawImage[y * imageW + x];
                    }
                    buff_Org[y - sy] /= (double)rc.Width;
                }
            }
            else if (nDir == DIR_EXFALL || nDir == DIR_EXRISE)
            {
                for (int y = ey - 1, nIndex = 0; y >= sy - 2; y--)
                {
                    for (int x = sx; x < ex; x++)
                    {
                        buff_Org[nIndex] += rawImage[y * imageW + x];
                    }
                    buff_Org[nIndex++] /= (double)rc.Width;
                }
            }

            return HC_EDGE_Get2ndDerivativeArrayFromLineBuff(buff_Org);
        }
        public static double[]/***/HC_EDGE_GetProjection_BTM_For_Hor_Derivative(byte[] rawImage, int imageW, int imageH, RectangleF rc, int nDir)
        {
            const int DIR_INFALL = 0;      // Direction = To Inside Falling
            const int DIR_INRISE = 1;      // Direction = to Inside Rising
            const int DIR_EXFALL = 2;      // Direction = to outside Falling
            const int DIR_EXRISE = 3;      // direction = to Outside Rising

            int sx = (int)rc.X;
            int ex = (int)rc.Width + sx;
            int sy = (int)rc.Y;
            int ey = (int)rc.Height + sy;

            double[] buff_Org = new double[(int)rc.Height + 2];

            if (nDir == DIR_INFALL || nDir == DIR_INRISE)
            {
                for (int y = ey - 1, nIndex = 0; y >= sy - 2; y--)
                {
                    for (int x = sx; x < ex; x++)
                    {
                        buff_Org[nIndex] += rawImage[y * imageW + x];
                    }
                    buff_Org[nIndex++] /= (double)rc.Width;
                }
            }
            else if (nDir == DIR_EXFALL || nDir == DIR_EXRISE)
            {
                for (int y = sy; y < ey + 2; y++)
                {
                    for (int x = sx; x < ex; x++)
                    {
                        buff_Org[y - sy] += rawImage[y * imageW + x];
                    }
                    buff_Org[y - sy] /= (double)rc.Width;
                }
            }

            return HC_EDGE_Get2ndDerivativeArrayFromLineBuff(buff_Org);
        }

        public static double[]/***/HC_EDGE_GetProjection_LFT_For_VER_Derivative(byte[] rawImage, int imageW, int imageH, RectangleF rc, int nDir)
        {
            const int DIR_INFALL = 0;      // Direction = To Inside Falling
            const int DIR_INRISE = 1;      // Direction = to Inside Rising
            const int DIR_EXFALL = 2;      // Direction = to outside Falling
            const int DIR_EXRISE = 3;      // direction = to Outside Rising

            int sx = (int)rc.X;
            int ex = (int)rc.Width + sx;
            int sy = (int)rc.Y;
            int ey = (int)rc.Height + sy;

            double[] buff_Org = new double[(int)rc.Width + 2];

            if (nDir == DIR_INFALL || nDir == DIR_INRISE)
            {
                for (int x = sx; x < ex + 2; x++)
                {
                    for (int y = sy; y < ey; y++)
                    {
                        buff_Org[x - sx] += rawImage[y * imageW + x];
                    }
                    buff_Org[x - sx] /= (double)rc.Width;
                }
            }
            else if (nDir == DIR_EXFALL || nDir == DIR_EXRISE)
            {
                for (int x = ex - 1, idx = 0; x >= sx - 2; x--)
                {
                    for (int y = sy; y < ey; y++)
                    {
                        buff_Org[idx] += rawImage[y * imageW + x];
                    }
                    buff_Org[idx++] /= (double)rc.Width;
                }
            }

            return HC_EDGE_Get2ndDerivativeArrayFromLineBuff(buff_Org);
        }
        public static double[]/***/HC_EDGE_GetProjection_RHT_For_VER_Derivative(byte[] rawImage, int imageW, int imageH, RectangleF rc, int nDir)
        {
            const int DIR_INFALL = 0;      // Direction = To Inside Falling
            const int DIR_INRISE = 1;      // Direction = to Inside Rising
            const int DIR_EXFALL = 2;      // Direction = to outside Falling
            const int DIR_EXRISE = 3;      // direction = to Outside Rising

            int sx = (int)rc.X;
            int ex = (int)rc.Width + sx;
            int sy = (int)rc.Y;
            int ey = (int)rc.Height + sy;

            double[] buff_Org = new double[(int)rc.Width + 2];

            if (nDir == DIR_INFALL || nDir == DIR_INRISE)
            {
                for (int x = ex - 1, idx = 0; x >= sx - 2; x--)
                {
                    for (int y = sy; y < ey; y++)
                    {
                        buff_Org[idx] += rawImage[y * imageW + x];
                    }
                    buff_Org[idx++] /= (double)rc.Width;
                }
            }
            else if (nDir == DIR_EXFALL || nDir == DIR_EXRISE)
            {
                for (int x = sx; x < ex + 2; x++)
                {
                    for (int y = sy; y < ey; y++)
                    {
                        buff_Org[x - sx] += rawImage[y * imageW + x];
                    }
                    buff_Org[x - sx] /= (double)rc.Width;
                }
            }

            return HC_EDGE_Get2ndDerivativeArrayFromLineBuff(buff_Org);
        }

        // 주어진 profile 값에서 min/max pos 결정 + 서브 픽셀링 하여, rectangle 값에 대입하여 방향에 따라 라인 결정해서 축에 대해 일괄 적용 후 반환
        public static List<PointF> HC_EDGE_ReplacePointList_Derivative_HOR_Average(double[] buff_Line, RectangleF rc, bool bTarget_TOP, int nDir)
        {
            const int DIR_INFALL = 0;      // Direction = To Inside Falling
            const int DIR_INRISE = 1;      // Direction = to Inside Rising
            const int DIR_EXFALL = 2;      // Direction = to outside Falling
            const int DIR_EXRISE = 3;      // direction = to Outside Rising

            List<PointF> listFiltered = new List<PointF>();

            int posMax = HC_ARRAY_GetMinElementPosition(buff_Line);
            int posMin = HC_ARRAY_GetMaxElementPosition(buff_Line);

            // get joint points positions
            int sx = (int)rc.X;
            int ex = (int)rc.Width + sx;
            int sy = (int)rc.Y;
            int ey = (int)rc.Height + sy;


            double fSubPixel = 0.0;

            if (bTarget_TOP == true)
            {
                if (nDir == DIR_INFALL)
                {
                    fSubPixel = HC_EDGE_GetSubPixelFromLineBuff(buff_Line, posMin);
                    for (int x = sx; x < ex; x++) { listFiltered.Add(new PointF(x, (float)(sy + posMin + fSubPixel))); }
                }
                else if (nDir == DIR_INRISE)
                {
                    fSubPixel = HC_EDGE_GetSubPixelFromLineBuff(buff_Line, posMax);
                    for (int x = sx; x < ex; x++) { listFiltered.Add(new PointF(x, (float)(sy + posMax + fSubPixel))); }
                }
                else if (nDir == DIR_EXRISE)
                {
                    fSubPixel = HC_EDGE_GetSubPixelFromLineBuff(buff_Line, posMax);
                    for (int x = sx; x < ex; x++) { listFiltered.Add(new PointF(x, (float)(ey - posMax - fSubPixel))); }
                }
                else if (nDir == DIR_EXFALL)
                {
                    fSubPixel = HC_EDGE_GetSubPixelFromLineBuff(buff_Line, posMin);
                    for (int x = sx; x < ex; x++) { listFiltered.Add(new PointF(x, (float)(ey - posMin - fSubPixel))); }
                }
            }
            else if (bTarget_TOP == false)
            {
                if (nDir == DIR_INRISE)
                {

                    fSubPixel = HC_EDGE_GetSubPixelFromLineBuff(buff_Line, posMax);
                    for (int x = sx; x < ex; x++) { listFiltered.Add(new PointF(x, (float)(ey - posMax - fSubPixel))); }
                }
                else if (nDir == DIR_INFALL)
                {
                    fSubPixel = HC_EDGE_GetSubPixelFromLineBuff(buff_Line, posMin);
                    for (int x = sx; x < ex; x++) { listFiltered.Add(new PointF(x, (float)(ey - posMin - fSubPixel))); }
                }
                else if (nDir == DIR_EXRISE)
                {
                    fSubPixel = HC_EDGE_GetSubPixelFromLineBuff(buff_Line, posMax);
                    for (int x = sx; x < ex; x++) { listFiltered.Add(new PointF(x, (float)(sy + posMax + fSubPixel))); }
                }
                else if (nDir == DIR_EXFALL)
                {
                    fSubPixel = HC_EDGE_GetSubPixelFromLineBuff(buff_Line, posMin);
                    for (int x = sx; x < ex; x++) { listFiltered.Add(new PointF(x, (float)(sy + posMin + fSubPixel))); }
                }
            }

            return listFiltered;
        }
        public static List<PointF> HC_EDGE_ReplacePointList_Derivative_VER_Average(double[] buff_Line, RectangleF rc, bool bTarget_LFT, int nDir)
        {
            const int DIR_INFALL = 0;      // Direction = To Inside Falling
            const int DIR_INRISE = 1;      // Direction = to Inside Rising
            const int DIR_EXFALL = 2;      // Direction = to outside Falling
            const int DIR_EXRISE = 3;      // direction = to Outside Rising

            List<PointF> listFiltered = new List<PointF>();

            int posMax = HC_ARRAY_GetMinElementPosition(buff_Line);
            int posMin = HC_ARRAY_GetMaxElementPosition(buff_Line);

            // get joint points positions
            int sx = (int)rc.X;
            int ex = (int)rc.Width + sx;
            int sy = (int)rc.Y;
            int ey = (int)rc.Height + sy;

            double fSubPixel = 0.0;

            if (bTarget_LFT == true)
            {
                if (nDir == DIR_INFALL)
                {
                    fSubPixel = HC_EDGE_GetSubPixelFromLineBuff(buff_Line, posMin);
                    for (int y = sy; y < ey; y++) { listFiltered.Add(new PointF((float)(sx + posMin + fSubPixel), y)); }
                }
                else if (nDir == DIR_INRISE)
                {
                    fSubPixel = HC_EDGE_GetSubPixelFromLineBuff(buff_Line, posMax);
                    for (int y = sy; y < ey; y++) { listFiltered.Add(new PointF((float)(sx + posMax + fSubPixel), y)); }
                }
                else if (nDir == DIR_EXFALL)
                {
                    fSubPixel = HC_EDGE_GetSubPixelFromLineBuff(buff_Line, posMin);
                    for (int y = sy; y < ey; y++) { listFiltered.Add(new PointF((float)(ex - posMin - fSubPixel), y)); }
                }
                else if (nDir == DIR_EXRISE)
                {
                    fSubPixel = HC_EDGE_GetSubPixelFromLineBuff(buff_Line, posMax);
                    for (int y = sy; y < ey; y++) { listFiltered.Add(new PointF((float)(ex - posMax - fSubPixel), y)); }
                }
            }
            else if (bTarget_LFT == false)
            {
                if (nDir == DIR_INFALL)
                {
                    fSubPixel = HC_EDGE_GetSubPixelFromLineBuff(buff_Line, posMin);
                    for (int y = sy; y < ey; y++) { listFiltered.Add(new PointF((float)(ex - posMin - fSubPixel), y)); }
                }
                else if (nDir == DIR_INRISE)
                {
                    fSubPixel = HC_EDGE_GetSubPixelFromLineBuff(buff_Line, posMax);
                    for (int y = sy; y < ey; y++) { listFiltered.Add(new PointF((float)(ex - posMax - fSubPixel), y)); }
                }
                else if (nDir == DIR_EXFALL)
                {
                    fSubPixel = HC_EDGE_GetSubPixelFromLineBuff(buff_Line, posMin);
                    for (int y = sy; y < ey; y++) { listFiltered.Add(new PointF((float)(sx + posMin + fSubPixel), y)); }
                }
                else if (nDir == DIR_EXRISE)
                {
                    fSubPixel = HC_EDGE_GetSubPixelFromLineBuff(buff_Line, posMax);
                    for (int y = sy; y < ey; y++) { listFiltered.Add(new PointF((float)(sx + posMax + fSubPixel), y)); }
                }
            }

            return listFiltered;
        }

        // get the every raw points based on directional 2nd derivative 170419 
        // this functions prepared for the point filtering or fitting which has outliers.
        public static List<PointF> HC_EDGE_GetPointList_Derivative_HOR(byte[] rawImage, int imageW, int imageH, RectangleF rc, bool bTarget_TOP, int nDir)
        {
            const int DIR_INFALL = 0;      // Direction = To Inside Falling
            const int DIR_INRISE = 1;      // Direction = to Inside Rising
            const int DIR_EXFALL = 2;      // Direction = to outside Falling
            const int DIR_EXRISE = 3;      // direction = to Outside Rising

            List<PointF> list = new List<PointF>();

            // get joint points positions
            int sx = (int)rc.X;
            int ex = (int)rc.Width + sx;
            int sy = (int)rc.Y;
            int ey = (int)rc.Height + sy;

            double fSubPixel = 0.0;

            double[] buff_Org = new double[(int)rc.Height + 2];

            if (bTarget_TOP == true)
            {
                #region FOR TOP REGION : IN-RISE & IN-FALL
                if (nDir == DIR_INFALL || nDir == DIR_INRISE)
                {
                    for (int x = sx; x < ex; x++)
                    {
                        Array.Clear(buff_Org, 0, buff_Org.Length);
                        for (int y = sy, nIndex = 0; y < ey + 2; y++)
                        {
                            buff_Org[nIndex++] = rawImage[y * imageW + x];
                        }
                        double[] buff_2nd = HC_EDGE_Get2ndDerivativeArrayFromLineBuff(buff_Org);

                        if (nDir == DIR_INFALL)
                        {
                            int maxPos = HC_ARRAY_GetMaxElementPosition(buff_2nd);
                            fSubPixel = HC_EDGE_GetSubPixelFromLineBuff(buff_2nd, maxPos);

                            list.Add(new PointF(x, (float)(sy + maxPos + fSubPixel)));
                        }
                        else if (nDir == DIR_INRISE)
                        {
                            int minPos = HC_ARRAY_GetMinElementPosition(buff_2nd);
                            fSubPixel = HC_EDGE_GetSubPixelFromLineBuff(buff_2nd, minPos);

                            list.Add(new PointF(x, (float)(sy + minPos + fSubPixel)));
                        }
                    }
                }
                #endregion

                #region FOR TOP REGION : IN-RISE AND IN-FALL

                else if (nDir == DIR_EXFALL || nDir == DIR_EXRISE)
                {
                    for (int x = sx; x < ex; x++)
                    {
                        Array.Clear(buff_Org, 0, buff_Org.Length);
                        for (int y = ey - 1, nIndex = 0; y >= sy - 2; y--)
                        {
                            buff_Org[nIndex++] = rawImage[y * imageW + x];
                        }
                        double[] buff_Top_2nd = HC_EDGE_Get2ndDerivativeArrayFromLineBuff(buff_Org);

                        if (nDir == DIR_EXFALL)
                        {
                            int maxPos = HC_ARRAY_GetMaxElementPosition(buff_Top_2nd);
                            fSubPixel = HC_EDGE_GetSubPixelFromLineBuff(buff_Top_2nd, maxPos);

                            list.Add(new PointF(x, (float)(ey - maxPos - fSubPixel)));
                        }
                        else if (nDir == DIR_EXRISE)
                        {
                            int minPos = HC_ARRAY_GetMinElementPosition(buff_Top_2nd);
                            fSubPixel = HC_EDGE_GetSubPixelFromLineBuff(buff_Top_2nd, minPos);
                            list.Add(new PointF(x, (float)(ey - minPos - fSubPixel)));
                        }
                    }
                }
                #endregion

            }
            else if (bTarget_TOP == false)
            {
                #region FOR BTM REGION : IN-RISE AND IN-FALL
                if (nDir == DIR_INFALL || nDir == DIR_INRISE)
                {
                    for (int x = sx; x < ex; x++)
                    {
                        Array.Clear(buff_Org, 0, buff_Org.Length);
                        for (int y = ey - 1, nIndex = 0; y >= sy - 2; y--)
                        {
                            buff_Org[nIndex++] = rawImage[y * imageW + x];
                        }
                        double[] buff_2nd = HC_EDGE_Get2ndDerivativeArrayFromLineBuff(buff_Org);

                        if (nDir == DIR_INFALL)
                        {
                            int minPos = HC_ARRAY_GetMaxElementPosition(buff_2nd);
                            fSubPixel = HC_EDGE_GetSubPixelFromLineBuff(buff_2nd, minPos);

                            list.Add(new PointF(x, (float)(ey - minPos - fSubPixel)));

                        }
                        else if (nDir == DIR_INRISE)
                        {
                            int maxPos = HC_ARRAY_GetMinElementPosition(buff_2nd);
                            fSubPixel = HC_EDGE_GetSubPixelFromLineBuff(buff_2nd, maxPos);

                            list.Add(new PointF(x, (float)(ey - maxPos - fSubPixel)));
                        }
                    }
                }
                #endregion

                #region FOR BTM REGION : EX-RISE AND EX-FALL
                else if (nDir == DIR_EXFALL || nDir == DIR_EXRISE)
                {
                    for (int x = sx; x < ex; x++)
                    {
                        Array.Clear(buff_Org, 0, buff_Org.Length);
                        for (int y = sy, nIndex = 0; y < ey + 2; y++)
                        {
                            buff_Org[nIndex++] = rawImage[y * imageW + x];
                        }
                        double[] buff_2nd = HC_EDGE_Get2ndDerivativeArrayFromLineBuff(buff_Org);

                        if (nDir == DIR_EXFALL)
                        {
                            int minPos = HC_ARRAY_GetMaxElementPosition(buff_2nd);
                            fSubPixel = HC_EDGE_GetSubPixelFromLineBuff(buff_2nd, minPos);
                            list.Add(new PointF(x, (float)(sy + minPos + fSubPixel)));

                        }
                        else if (nDir == DIR_EXRISE)
                        {
                            int maxPos = HC_ARRAY_GetMinElementPosition(buff_2nd);
                            fSubPixel = HC_EDGE_GetSubPixelFromLineBuff(buff_2nd, maxPos);
                            list.Add(new PointF(x, (float)(sy + maxPos + fSubPixel)));
                        }
                    }
                }
                #endregion
            }

            return list;
        }
        public static List<PointF> HC_EDGE_GetPointList_Derivative_VER(byte[] rawImage, int imageW, int imageH, RectangleF rc, bool bTarget_LFT, int nDir)
        {

            const int DIR_INFALL = 0;      // Direction = To Inside Falling
            const int DIR_INRISE = 1;      // Direction = to Inside Rising
            const int DIR_EXFALL = 2;      // Direction = to outside Falling
            const int DIR_EXRISE = 3;      // direction = to Outside Rising

            List<PointF> list = new List<PointF>();

            // get joint points positions
            int sx = (int)rc.X;
            int ex = (int)rc.Width + sx;
            int sy = (int)rc.Y;
            int ey = (int)rc.Height + sy;

            double fSubPixel = 0.0;

            double[] buff_Org = new double[(int)rc.Width + 2];

            if (bTarget_LFT == true)
            {
                #region FOR TOP REGION : IN-RISE & IN-FALL
                if (nDir == DIR_INFALL || nDir == DIR_INRISE)
                {
                    for (int y = sy; y < ey; y++)
                    {
                        Array.Clear(buff_Org, 0, buff_Org.Length);
                        for (int x = sx, nIndex = 0; x < ex + 2; x++)
                        {
                            buff_Org[nIndex++] = rawImage[y * imageW + x];
                        }
                        double[] buff_2nd = HC_EDGE_Get2ndDerivativeArrayFromLineBuff(buff_Org);

                        if (nDir == DIR_INFALL)
                        {
                            int maxPos = HC_ARRAY_GetMaxElementPosition(buff_2nd);
                            fSubPixel = HC_EDGE_GetSubPixelFromLineBuff(buff_2nd, maxPos);
                            list.Add(new PointF((float)(sx + maxPos + fSubPixel), y));
                        }
                        else if (nDir == DIR_INRISE)
                        {
                            int minPos = HC_ARRAY_GetMinElementPosition(buff_2nd);
                            fSubPixel = HC_EDGE_GetSubPixelFromLineBuff(buff_2nd, minPos);
                            list.Add(new PointF((float)(sx + minPos + fSubPixel), y));
                        }

                    }
                }
                #endregion

                #region FOR TOP REGION : IN-RISE AND IN-FALL
                else if (nDir == DIR_EXFALL || nDir == DIR_EXRISE)
                {
                    for (int y = sy; y < ey; y++)
                    {
                        Array.Clear(buff_Org, 0, buff_Org.Length);
                        for (int x = ex - 1, nIndex = 0; x >= sx - 2; x--)
                        {
                            buff_Org[nIndex++] = rawImage[y * imageW + x];
                        }
                        double[] buff_2nd = HC_EDGE_Get2ndDerivativeArrayFromLineBuff(buff_Org);

                        if (nDir == DIR_EXFALL)
                        {
                            int maxPos = HC_ARRAY_GetMaxElementPosition(buff_2nd);
                            fSubPixel = HC_EDGE_GetSubPixelFromLineBuff(buff_2nd, maxPos);
                            list.Add(new PointF((float)(ex - maxPos - fSubPixel), y));
                        }
                        else if (nDir == DIR_EXRISE)
                        {
                            int minPos = HC_ARRAY_GetMinElementPosition(buff_2nd);
                            fSubPixel = HC_EDGE_GetSubPixelFromLineBuff(buff_2nd, minPos);
                            list.Add(new PointF((float)(ex - minPos - fSubPixel), y));
                        }
                    }
                }
                #endregion
            }
            else if (bTarget_LFT == false)
            {
                #region FOR BTM REGION : IN-RISE AND IN-FALL

                if (nDir == DIR_INFALL || nDir == DIR_INRISE)
                {
                    for (int y = sy; y < ey; y++)
                    {
                        Array.Clear(buff_Org, 0, buff_Org.Length);
                        for (int x = ex - 1, nIndex = 0; x >= sx - 2; x--)
                        {
                            buff_Org[nIndex++] = rawImage[y * imageW + x];
                        }
                        double[] buff_RHT_2nd = HC_EDGE_Get2ndDerivativeArrayFromLineBuff(buff_Org);

                        if (nDir == DIR_INFALL)
                        {

                            int maxPos = HC_ARRAY_GetMaxElementPosition(buff_RHT_2nd);
                            fSubPixel = HC_EDGE_GetSubPixelFromLineBuff(buff_RHT_2nd, maxPos);
                            list.Add(new PointF((float)(ex - maxPos - fSubPixel), y));
                        }
                        else if (nDir == DIR_INRISE)
                        {
                            int minPos = HC_ARRAY_GetMinElementPosition(buff_RHT_2nd);
                            fSubPixel = HC_EDGE_GetSubPixelFromLineBuff(buff_RHT_2nd, minPos);
                            list.Add(new PointF((float)(ex - minPos - fSubPixel), y));
                        }
                    }
                }
                #endregion

                #region FOR BTM REGION : EX-RISE AND EX-FALL
                else if (nDir == DIR_EXFALL || nDir == DIR_EXRISE)
                {
                    for (int y = sy; y < ey; y++)
                    {
                        Array.Clear(buff_Org, 0, buff_Org.Length);
                        for (int x = sx, nIndex = 0; x < ex + 2; x++)
                        {
                            buff_Org[nIndex++] = rawImage[y * imageW + x];
                        }
                        double[] buff_2nd = HC_EDGE_Get2ndDerivativeArrayFromLineBuff(buff_Org);
                        if (nDir == DIR_EXFALL)
                        {
                            int maxPos = HC_ARRAY_GetMaxElementPosition(buff_2nd);
                            fSubPixel = HC_EDGE_GetSubPixelFromLineBuff(buff_2nd, maxPos);
                            list.Add(new PointF((float)(sx + maxPos + fSubPixel), y));
                        }
                        else if (nDir == DIR_EXRISE)
                        {
                            int minPos = HC_ARRAY_GetMinElementPosition(buff_2nd);
                            fSubPixel = HC_EDGE_GetSubPixelFromLineBuff(buff_2nd, minPos);

                            list.Add(new PointF((float)(sx + minPos + fSubPixel), y));
                        }
                    }
                }
                #endregion
            }
            return list;
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

        //*******************************************************************************************
        // Zerocrossing based Edge detection
        //*******************************************************************************************

        /// <summary> 170426 : cd-meter project
        /// GetThe zerocrossing point from the line buff
        /// </summary>
        /// <param name="buff"></param>
        /// <param name="sign"></param>
        /// <param name="mag"></param>
        /// <returns></returns>
        public static double HC_EDGE_GetZeroCrossingPointFromLineBuff(byte [] buff, int sign, double mag )
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

        // for overlay 방향에 따른 제로크로싱 리스트를 계산
        public static List<PointF> HC_EDGE_GetRawPoints_Hor_ZeroCross(byte[] rawImage, int imageW, int imageH, RectangleF rc, bool bTarget_TOP, int nDir, double mag)
        {
            // get joint points positions
            int sx = (int)rc.X;
            int ex = (int)rc.Width + sx;
            int sy = (int)rc.Y;
            int ey = (int)rc.Height + sy;

            byte[] buffLine = new byte[ey - sy];

            List<PointF> list = new List<PointF>();

            // top side reverse 
            if (bTarget_TOP == true)
            {
                for (int x = sx; x < ex; x++)
                {
                    Array.Clear(buffLine, 0, buffLine.Length);
                    for (int y = ey - 1, nIndex = 0; y >= sy; y--)
                    {
                        buffLine[nIndex++] = rawImage[y * imageW + x];
                    }

                    double fSubPos = HC_EDGE_GetZeroCrossingPointFromLineBuff(buffLine, nDir, mag);
                    list.Add(new PointF(x, ey - (float)fSubPos));
                }
            }
            // btm side 
            else if (bTarget_TOP == false)
            {
                for (int x = sx; x < ex; x++)
                {
                    Array.Clear(buffLine, 0, buffLine.Length);
                    for (int y = sy, nIndex = 0; y < ey; y++)
                    {
                        buffLine[nIndex++] = rawImage[y * imageW + x];
                    }
                    double fSubPos = HC_EDGE_GetZeroCrossingPointFromLineBuff(buffLine, nDir, mag);
                    list.Add(new PointF(x, sy + (float)fSubPos));
                }
            }
            return list;
        }
        public static List<PointF> HC_EDGE_GetRawPoints_VER_ZeroCross(byte[] rawImage, int imageW, int imageH, RectangleF rc, bool bTarget_LFT, int nDir, double mag)
        {
            // get joint points positions
            int sx = (int)rc.X;
            int ex = (int)rc.Width + sx;
            int sy = (int)rc.Y;
            int ey = (int)rc.Height + sy;

            byte[] buffLine = new byte[ex - sx];

            List<PointF> list = new List<PointF>();

            // top side reverse 
            if (bTarget_LFT == true)
            {
                for (int y = sy; y < ey; y++)
                {
                    Array.Clear(buffLine, 0, buffLine.Length);
                    for (int x = ex - 1, nIndex = 0; x >= sx; x--)
                    {
                        buffLine[nIndex++] = rawImage[y * imageW + x];
                    }
                    double fSubPos = HC_EDGE_GetZeroCrossingPointFromLineBuff(buffLine, nDir, mag);
                    list.Add(new PointF(ex - (float)fSubPos, y));
                }
            }
            // btm side 
            else if (bTarget_LFT == false)
            {
                for (int y = sy; y < ey; y++)
                {
                    Array.Clear(buffLine, 0, buffLine.Length);
                    for (int x = sx, nIndex = 0; x < ex; x++)
                    {
                        buffLine[nIndex++] = rawImage[y * imageW + x];
                    }

                    double fSubPos = HC_EDGE_GetZeroCrossingPointFromLineBuff(buffLine, nDir, mag);
                    list.Add(new PointF(sx + (float)fSubPos, y));
                }
            }
            return list;
        }


        //*******************************************************************************************
        // Laplacian of Gaussian based Edge detection
        //*******************************************************************************************
        /// <summary> 170426
        /// 라인버퍼로부터 Min/Max position이 결정된경우 라인버퍼에서 해당 위치에 대해서 전후 subpixeling 하기 : from the LOG Edge detection 
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="lineBuff"></param>
        /// <param name="nPos"></param>
        /// <returns></returns>
        public static double HC_EDGE_GetSubPixelFromLineBuff<T>(T[] lineBuff, int nPos)
        {
            double pa = 0; double pb = 0; double pc = 0;

            double fSubPixel = 0;

            try
            {
                if (nPos != 0 && nPos < lineBuff.Length - 1)
                {
                    pa = Convert.ToDouble(lineBuff[nPos - 1]);
                    pb = Convert.ToDouble(lineBuff[nPos + 0]);
                    pc = Convert.ToDouble(lineBuff[nPos + 1]);

                    // simple quadratic interpolation
                    fSubPixel = 0.5 * (pa - pc) / (pa - (2 * pb) + pc);
                }
            }
            catch (Exception ex)
            {
                Console.WriteLine(ex.ToString());
            }


            return fSubPixel;
        }
        public static double HC_EDGE_GetLogPos(byte[] rawImage, int imageW, int imageH, PointF[] arrPoints, int nSign)
        {
            double[] fKernel = HC_FILTER_GetLogKernel(9);

            double[] fImage = new double[arrPoints.Length];

            int KSIZE = (int)Math.Sqrt(fKernel.Length);
            int GAP = KSIZE / 2;

            for (int i = 0; i < arrPoints.Length; i++)
            {
                int x = (int)arrPoints.ElementAt(i).X;
                int y = (int)arrPoints.ElementAt(i).Y;

                double kernelSum = 0;
                for (int j = -GAP; j <= GAP; j++)
                {
                    for (int k = -GAP; k <= GAP; k++)
                    {
                        kernelSum += (fKernel[(j + GAP) * KSIZE + k + GAP] * rawImage[(y - j) * imageW + (x - k)]);
                    }
                }
                fImage[i] = kernelSum;
            }

            double fValue = 0;
            int nPos = 0;

            /***/
            if (nSign == -1) { fValue = fImage.Min(); }
            else if (nSign == +1) { fValue = fImage.Max(); }
            nPos = Array.IndexOf(fImage, fValue);

            double fSubPos = HC_EDGE_GetSubPixelFromLineBuff(fImage, nPos);

            return nPos + fSubPos;

        }
        public static double HC_EDGE_GetLoG_PosMax(byte[] rawImage, int imageW, int imageH, PointF[] list)
        {
            double[] arrDerivative = new double[list.Length];

            for (int i = 0; i < list.Length; i++)
            {
                int xx = (int)list.ElementAt(i).X;
                int yy = (int)list.ElementAt(i).Y;

                arrDerivative[i] = rawImage[yy * imageW + xx + 1] + rawImage[yy * imageW + xx - 1] - (2 * rawImage[yy * imageW + xx]);
            }

            double fMaxVal = arrDerivative.Max();
            int nMaxPos = Array.IndexOf(arrDerivative, fMaxVal);

            double fSubPixel = HC_EDGE_GetSubPixelFromLineBuff(arrDerivative, nMaxPos);

            return nMaxPos + fSubPixel;
        }
        public static double HC_EDGE_GetLoG_PosMin(byte[] rawImage, int imageW, int imageH, PointF[] list)
        {
            double[] arrDerivative = new double[list.Length];

            for (int i = 0; i < list.Length; i++)
            {
                int xx = (int)list.ElementAt(i).X;
                int yy = (int)list.ElementAt(i).Y;

                arrDerivative[i] = rawImage[(yy + 0) * imageW + xx + 1] +
                                   rawImage[(yy + 0) * imageW + xx - 1] +
                                   rawImage[(yy + 1) * imageW + xx + 0] +
                                   rawImage[(yy - 1) * imageW + xx + 0] -
                                   (4 * rawImage[yy * imageW + xx]);
            }

            double fMin = arrDerivative.Min();
            int nMinPos = Array.IndexOf(arrDerivative, fMin);

            double fSubpixel = HC_EDGE_GetSubPixelFromLineBuff(arrDerivative, nMinPos);

            return nMinPos + fSubpixel;
        }

        public static List<PointF> GetRawPoints_HOR_LOG(byte[] rawImage, int imageW, int imageH, RectangleF rc, bool bTarget_TOP, int nDir)
        {
            // get joint points positions
            int sx = (int)rc.X;
            int ex = (int)rc.Width + sx;
            int sy = (int)rc.Y;
            int ey = (int)rc.Height + sy;

            PointF[] buffPoints = new PointF[(int)rc.Height];

            List<PointF> list = new List<PointF>();

            // top side reverse 
            if (bTarget_TOP == true)
            {
                for (int x = sx; x < ex; x++)
                {
                    Array.Clear(buffPoints, 0, buffPoints.Length);
                    for (int y = ey - 1, nIndex = 0; y >= sy; y--)
                    {
                        buffPoints[nIndex++] = new PointF(x, y);
                    }

                    double fSubPos = HC_EDGE_GetLogPos(rawImage, imageW, imageH, buffPoints, nDir);
                    list.Add(new PointF(x, ey - (float)fSubPos));
                }
            }
            // btm side 
            else if (bTarget_TOP == false)
            {
                for (int x = sx; x < ex; x++)
                {
                    Array.Clear(buffPoints, 0, buffPoints.Length);
                    for (int y = sy, nIndex = 0; y < ey; y++)
                    {
                        buffPoints[nIndex++] = new PointF(x, y);
                    }

                    double fSubPos = HC_EDGE_GetLogPos(rawImage, imageW, imageH, buffPoints, nDir);
                    list.Add(new PointF(x, sy + (float)fSubPos));
                }

            }
            return list;
        }
        public static List<PointF> GetRawPoints_VER_LOG(byte[] rawImage, int imageW, int imageH, RectangleF rc, bool bTarget_LFT, int nDir)
        {
            // get joint points positions
            int sx = (int)rc.X;
            int ex = (int)rc.Width + sx;
            int sy = (int)rc.Y;
            int ey = (int)rc.Height + sy;

            PointF[] buffPoints = new PointF[(int)rc.Width];

            List<PointF> list = new List<PointF>();

            // top side reverse 
            if (bTarget_LFT == true)
            {
                for (int y = sy; y < ey; y++)
                {
                    Array.Clear(buffPoints, 0, buffPoints.Length);
                    for (int x = ex - 1, nIndex = 0; x >= sx; x--)
                    {
                        buffPoints[nIndex++] = new PointF(x, y);
                    }
                    double fSubPos = HC_EDGE_GetLogPos(rawImage, imageW, imageH, buffPoints, nDir);

                    list.Add(new PointF(ex - (float)fSubPos, y));
                }
            }
            // btm side 
            else if (bTarget_LFT == false)
            {
                for (int y = sy; y < ey; y++)
                {
                    Array.Clear(buffPoints, 0, buffPoints.Length);
                    for (int x = sx, nIndex = 0; x < ex; x++)
                    {
                        buffPoints[nIndex++] = new PointF(x, y);
                    }
                    double fSubPos = HC_EDGE_GetLogPos(rawImage, imageW, imageH, buffPoints, nDir);
                    list.Add(new PointF(sx + (float)fSubPos, y));
                }
            }
            return list;
        }
        /// <summary> 170426
        /// If there is floating position, make the bilinear interpolation value
        /// </summary>
        /// <param name="rawImage"></param>
        /// <param name="imageW"></param>
        /// <param name="imageH"></param>
        /// <param name="px"></param>
        /// <param name="py"></param>
        /// <returns></returns>
        public static double HC_GetInterpolatedValue(byte[] rawImage, int imageW, int imageH, float px, float py)
        {
            double cx = px;
            double cy = py;
            int x1 = (int)Math.Floor(cx);
            int x2 = (int)Math.Ceiling(cx);
            int y1 = (int)Math.Floor(cy);
            int y2 = (int)Math.Ceiling(cy);

            int q11 = rawImage[y1 * imageW + x1];
            int q12 = rawImage[y2 * imageW + x1];
            int q21 = rawImage[y1 * imageW + x2];
            int q22 = rawImage[y2 * imageW + x2];

            double fInterplated = _GetBilinearInterpolation(cx, cy, x1, x2, y1, y2, q11, q12, q21, q22);
            return fInterplated;
        }
        public static double _GetBilinearInterpolation(double cx, double cy, double x1, double x2, double y1, double y2, double q11, double q12, double q21, double q22)
        {
            double r1 = (((x2 - cx) / (x2 - x1)) * q11) + (((cx - x1) / (x2 - x1)) * q21);
            double r2 = (((x2 - cx) / (x2 - x1)) * q12) + (((cx - x1) / (x2 - x1)) * q22);
            double pvalue = (((y2 - cy) / (y2 - y1)) * r1) + (((cy - y1) / (y2 - y1)) * r2);
            return pvalue;
        }

        //*******************************************************************************************
        // Prewitt based Edge detection 
        //*******************************************************************************************

        // 단순 Prewitt 연산 array 만 받는다 
        public static double[] GetPrewitBuffLine(byte[] buffLine, int nMeasureType)
        {
            const int RISING = 0;
            const int FALLING = 1;

            double[] profile = new double[buffLine.Length];

            // accumulation for each position
            for (int nIndex = 1; nIndex < buffLine.Length - 1; nIndex++)
            {
                if (nMeasureType == RISING)
                {
                    profile[nIndex] += buffLine[(nIndex + 1)] - buffLine[(nIndex - 1)];
                }
                else if (nMeasureType == FALLING)
                {
                    profile[nIndex] += buffLine[(nIndex - 1)] - buffLine[(nIndex + 1)];
                }
            }

            return profile;
        }
        //  capsulated Edge detection function for only horizontal and vertical  170412 
        public static List<PointF> GetRawPoints_HOR_Prewitt(byte[] rawImage, int imageW, int imageH, RectangleF rc, bool bTarget_TOP, int nDir)
        {
            // get joint points positions
            int sx = (int)rc.X;
            int ex = (int)rc.Width + sx;
            int sy = (int)rc.Y;
            int ey = (int)rc.Height + sy;

            byte[] buffLine = new byte[ey - sy];

            List<PointF> list = new List<PointF>();

            // top side reverse 
            if (bTarget_TOP == true)
            {
                for (int x = sx; x < ex; x++)
                {
                    Array.Clear(buffLine, 0, buffLine.Length);
                    for (int y = sy, nIndex = 0; y < ey; y++)
                    {
                        buffLine[nIndex++] = rawImage[y * imageW + x];
                    }
                    double[] projH = _GetPrewitBuffLine(buffLine, nDir);
                    float yy = (float)_GetNewtonRapRes(projH);

                    list.Add(new PointF(x, (float)(sy + yy)));
                }
            }
            // btm side 
            else if (bTarget_TOP == false)
            {
                for (int x = sx; x < ex; x++)
                {
                    Array.Clear(buffLine, 0, buffLine.Length);
                    for (int y = ey - 1, nIndex = 0; y >= sy; y--)
                    {
                        buffLine[nIndex++] = rawImage[y * imageW + x];
                    }
                    double[] projH = GetPrewitBuffLine(buffLine, nDir);
                    float yy = (float)_GetNewtonRapRes(projH);

                    list.Add(new PointF(x, (float)(ey + yy)));
                }
            }
            return list;
        }
        public static List<PointF> GetRawPoints_VER_Prewitt(byte[] rawImage, int imageW, int imageH, RectangleF rc, bool bTarget_LFT, int nDir)
        {
            // get joint points positions
            int sx = (int)rc.X;
            int ex = (int)rc.Width + sx;
            int sy = (int)rc.Y;
            int ey = (int)rc.Height + sy;

            byte[] buffLine = new byte[ex - sx];

            List<PointF> list = new List<PointF>();

            // top side reverse 
            if (bTarget_LFT == true)
            {

                for (int y = sy; y < ey; y++)
                {
                    Array.Clear(buffLine, 0, buffLine.Length);
                    for (int x = ex - 1, nIndex = 0; x >= sx; x--)
                    {
                        buffLine[nIndex++] = rawImage[y * imageW + x];
                    }
                    double[] projV = GetPrewitBuffLine(buffLine, nDir);
                    float xx = (float)_GetNewtonRapRes(projV);

                    list.Add(new PointF(ex - xx, y));
                }
            }
            // btm side 
            else if (bTarget_LFT == false)
            {
                for (int y = sy; y < ey; y++)
                {
                    for (int x = sx, nIndex = 0; x < ex; x++)
                    {
                        buffLine[nIndex++] = rawImage[y * imageW + x];
                    }
                    double[] projV = GetPrewitBuffLine(buffLine, nDir);
                    float xx = (float)_GetNewtonRapRes(projV);

                    list.Add(new PointF(sx + xx, y));
                }
            }
            return list;
        }
        public static double[] _GetPrewitBuffLine(byte[] buffLine, int nMeasureType)
        {
            const int RISING = 0;
            const int FALLING = 1;

            double[] profile = new double[buffLine.Length];

            // accumulation for each position
            for (int nIndex = 1; nIndex < buffLine.Length - 1; nIndex++)
            {
                if (nMeasureType == RISING)
                {
                    profile[nIndex] += buffLine[(nIndex + 1)] - buffLine[(nIndex - 1)];
                }
                else if (nMeasureType == FALLING)
                {
                    profile[nIndex] += buffLine[(nIndex - 1)] - buffLine[(nIndex + 1)];
                }
            }

            return profile;
        }
        public static double _NewtonRaphson(double[] pNewton, int szPoly, double tStart)
        {
            double tSlope = 0;
            double tPoint = tStart;

            for (int itr = 0; itr < 10; itr++)
            {
                double tValue = tSlope = 0.0;

                for (int i = 0; i < szPoly; i++)
                {
                    tValue += pNewton[i] * Math.Pow(tPoint, (double)(szPoly - 1 - i));
                }
                for (int i = 0; i < szPoly - 1; i++)
                {
                    tSlope += (double)(szPoly - 1 - i) * pNewton[i] * Math.Pow(tPoint, (double)(szPoly - 2 - i));
                }
                double bPoint = tPoint;

                if (Math.Abs(tSlope) < 1e-10) break;
                if (Math.Abs(tValue) < 1e-16) break;

                tPoint = (tSlope * tPoint - tValue) / tSlope;
                if (Math.Abs(bPoint - tPoint) < 1e-5) break;
            }
            return tPoint;
        }
        public static float _GetNewtonRapRes(double[] arrAccProfile)
        {
            // set kernel size by default 
            const int KSIZE = 5;
            double[][] pMatrix = new double[5][];

            // allocation
            for (int y = 0; y < KSIZE; y++) { pMatrix[y] = new double[6]; }

            // get max position
            double fMin = arrAccProfile.Min();
            int nMinPos = Array.IndexOf(arrAccProfile, fMin);

            // assignment
            for (int y = 0; y < KSIZE; y++)
            {
                for (int x = 0; x < KSIZE; x++)
                {
                    pMatrix[y][x] = Math.Pow((double)(y + 1), (double)(KSIZE - 1 - x));
                }
                int nPos = (int)(nMinPos - ((KSIZE - 1) / 2.0) + y);

                // additive fuck exception : incase of Position want to place on the asshole 
                if (nPos > 0 && arrAccProfile.Length > nPos)
                {
                    pMatrix[y][KSIZE] = arrAccProfile[nPos];
                }
                else
                {
                    pMatrix[y][KSIZE] = 0;
                }
            }

            // fucking gauss
            _GaussElimination(pMatrix, KSIZE + 1, KSIZE);

            double[] pNewton = new double[KSIZE];

            for (int x = 0; x < KSIZE - 1; x++)
            {
                pNewton[x] = (double)(KSIZE - 1 - x) * pMatrix[x][KSIZE];
            }

            double NRValue = _NewtonRaphson(pNewton, KSIZE - 1, (double)(KSIZE + 1) / 2.0);

            // set fucking value
            return Convert.ToSingle(nMinPos + NRValue);
        }
        public static void _GaussElimination(double[][] pMatrix, int szX, int szY)
        {
            // X must be bigger than Y by 1

            // Left Diagonal
            for (int i = 0; i < szX - 2; i++)
            {
                for (int j = i + 1; j < szY; j++)
                {
                    if (pMatrix[j][i] * pMatrix[i][i] != double.NaN)
                    {
                        double eCoeff = pMatrix[j][i] / pMatrix[i][i];
                        for (int k = i; k < szX; k++)
                        {
                            pMatrix[j][k] = pMatrix[j][k] - pMatrix[i][k] * eCoeff;
                        }
                    }
                }
            }

            // Right Diagonal
            for (int j = 0; j < szY - 1; j++)
            {
                for (int i = j + 1; i < szX - 1; i++)
                {
                    if (pMatrix[j][i] * pMatrix[i][i] != double.NaN)
                    {
                        double eCoeff = pMatrix[j][i] / pMatrix[i][i];
                        for (int k = i; k < szX; k++)
                        {
                            pMatrix[j][k] = pMatrix[j][k] - pMatrix[i][k] * eCoeff;
                        }
                    }
                }
            }

            for (int j = 0; j < szY; j++)
            {
                if (pMatrix[j][j] != double.NaN)
                {
                    double eCoeff = pMatrix[j][j];
                    pMatrix[j][j] = 1.0;
                    pMatrix[j][szX - 1] = pMatrix[j][szX - 1] / eCoeff;
                }
            }
        }

        //********************************************************************************************************************************
        // Point Filtering
        //********************************************************************************************************************************

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
