using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Drawing;
using System.Drawing.Imaging;

using DispObject;

namespace ReinlessLib
{
    public static partial class CReinlessLib
    {
        public static byte[] HC_Fill_Polar_Raw(byte[] rawImage, int imageW, int imageH, int nRadius, int nValue)
        {
            int nCenterX = imageW / 2;
            int nCenterY = imageH / 2;

            byte[] rawRes = new byte[imageW * imageH];

            for (int y = 0; y < imageH; y++)
            {
                for (int x = 0; x < imageW; x++)
                {
                    if (Math.Sqrt(Math.Pow((nCenterX - x), 2.0) + Math.Pow((nCenterY - y), 2.0)) > nRadius)
                    {
                        rawRes[y * imageW + x] = Convert.ToByte(nValue);
                    }
                }
            }

            return rawRes;
        }
        public static byte[] HC_Fill_Henning(byte [] rawImage, int imageW, int imageH)
        {
            byte[] rawImage_Pre = new byte[imageW * imageH];

            double X, Y, Radius;
            double x0 = imageW / 2;
            double y0 = imageH / 2;

            //double Sigma = 130;
            double re; // distance from edge
            double x1, y1; // circle edge coord
            double b = 0.003; //decay factor. this defines edge-effectness
            double CE = 10; //How edge will be darken

            x1 = y1 = 0.95 * imageW / 2;

            double rt = 30;

            for (int y = 0; y < imageH; y++)
            {
                for (int x = 0; x < imageW; x++)
                {
                    X = x - imageW / 2;
                    Y = y - imageH / 2;
                    Radius = Math.Sqrt(X * X + Y * Y);

                    re = x1 - Radius;

                    if (Radius < (0.95 * imageW / 2))
                    {
                        if (re > rt) rawImage_Pre[y * imageH + x] = rawImage[y * imageH + x];
                        else rawImage_Pre[y * imageH + x] = (byte)Math.Round(CE + (rawImage[y * imageH + x] - CE) * (-1 + Math.Exp(-b * re)) / (-1 + Math.Exp(-b * rt)));

                        if (rawImage_Pre[y * imageH + x] < 0) rawImage_Pre[y * imageH + x] = 0;
                    }
                    else rawImage_Pre[y * imageH + x] = (byte)0;

                }
            }

            return rawImage_Pre;
        }

        public static void HC_FILL_Digonal_Region(byte[] rawImage, int imageW, int imageH, int nDirection, int nWidth, int nHeight, int nValue)
        {
            if (nDirection == 0)    // left top
            {
                for (int x = 0; x < nWidth; x++) { for (int y = 0; y < nHeight; y++) { rawImage[y * imageW + x] = (Byte)nValue; } }
            }
            else if (nDirection == 1) // right top
            {
                for (int x = imageW - nWidth; x < imageW; x++) { for (int y = 0; y < nHeight; y++) { rawImage[y * imageW + x] = (Byte)nValue; } }
            }
            else if (nDirection == 2) // left btm
            {
                for (int x = 0; x < nWidth; x++) { for (int y = imageH - nHeight; y < imageH; y++) { rawImage[y * imageW + x] = (Byte)nValue; } }
            }
            else if (nDirection == 3) // right btm
            {
                for (int x = imageW - nWidth; x < imageW; x++) { for (int y = imageH - nHeight; y < imageH; y++) { rawImage[y * imageW + x] = (Byte)nValue; } }
            }
        }

        public static byte[] HC_FILL_Boundary(byte[] rawImage, int imageW, int imageH, int nGap, int nValue)
        {
            // top side
            for (int y = 0; y < nGap; y++)
            {
                for (int x = 0; x < imageW; x++){rawImage[y * imageW + x] = (byte)nValue;}
            }

            // btm side 
            for (int y = imageH - 1; y >= imageH - nGap; y--)
            {
                for (int x = 0; x < imageW; x++){rawImage[y * imageW + x] = (byte)nValue;}
            }

            // left side
            for (int x = 0; x < nGap; x++)
            {
                for (int y = 0; y < imageH; y++){rawImage[y * imageW + x] = (byte)nValue;}
            }

            // right side
            for (int x = imageW - 1; x >= imageW - nGap; x--)
            {
                for (int y = 0; y < imageH; y++){rawImage[y * imageW + x] = (byte)nValue;}
            }
            return rawImage;
        }
        public static byte[] HC_FILL_Boundary(byte[] rawImage, int imageW, int imageH, double fRatio, byte value)
        {
            int nGap = 0;
            nGap = Convert.ToInt32(imageH * fRatio);

            // Boundary TOP Side
            for (int y = 0; y < nGap; y++) { for (int x = 0; x < imageW; x++) { rawImage[y * imageW + x] = value; } }

            nGap = Convert.ToInt32(imageH * (1 - fRatio));
            // Boundary Btm Side
            for (int y = nGap; y < imageH; y++) { for (int x = 0; x < imageW; x++) { rawImage[y * imageW + x] = value; } }

            nGap = Convert.ToInt32(imageW * fRatio);
            // Boundary Left Side
            for (int y = 0; y < imageH; y++) { for (int x = 0; x < nGap; x++) { rawImage[y * imageW + x] = value; } }

            nGap = Convert.ToInt32(imageW * (1 - fRatio));
            // Boundary Right Side
            for (int y = 0; y < imageH; y++) { for (int x = nGap; x < imageW; x++) { rawImage[y * imageW + x] = value; } }
            return rawImage;
        }
        public static byte[] HC_Fill_Discover_Bolder_Adjacent_RAW(byte[] rawImage, int nWidth, int nHeight, int nGap, int nValue)
        {
            // Boundary TOP Side
            int reverseGap = nGap * 2;
            for (int y = 0; y < nGap; y++) { for (int x = nGap; x < nWidth - nGap; x++) { rawImage[y * nWidth + x] = rawImage[(y + reverseGap) * nWidth + x]; } reverseGap--; }

            // Boundary Btm Side
            reverseGap = nGap * 2;
            for (int y = nHeight - nGap; y < nHeight; y++) { for (int x = nGap; x < nWidth - nGap; x++) { rawImage[y * nWidth + x] = rawImage[(y - reverseGap) * nWidth + x]; } reverseGap--; }

            // Boundary Left Side
            reverseGap = nGap * 2;
            for (int x = 0; x < nGap; x++) { for (int y = nGap; y < nHeight - nGap; y++) { rawImage[y * nWidth + x] = rawImage[y * nWidth + (x + reverseGap)]; } reverseGap--; }

            // Boundary Right Side
            reverseGap = nGap * 2;
            for (int x = 0; x < nGap; x++) { for (int y = nGap; y < nHeight - nGap; y++) { rawImage[y * nWidth + x] = rawImage[y * nWidth + (x - reverseGap)]; } reverseGap--; }
            return rawImage;
        }

        public static void HC_FILL_Region_Horizontal(byte[] rawImage, int imageW, int imageH, bool bLeft, bool bRight, double fCutRatio, byte value)
        {
            // position exception
            if (fCutRatio <= 0 || fCutRatio > 1.0) return;

            int nCutPos = 0;

            if (bLeft)
            {
                nCutPos = Convert.ToInt32(imageW * fCutRatio);

                Parallel.For(0, imageH, y =>
                {
                    // cut from 0 to cutPos
                    for (int x = 0; x < nCutPos; x++)
                    {
                        rawImage[y * imageW + x] = value;
                    }
                });
            }

            if (bRight)
            {
                nCutPos = Convert.ToInt32(imageW * (1 - fCutRatio));
                Parallel.For(0, imageH, y =>
                {
                    // cut from cutpos to endPos
                    for (int x = nCutPos; x < imageW; x++)
                    {
                        rawImage[y * imageW + x] = value;
                    }
                });
            }

        }
        public static void HC_FILL_Region_Vertical(byte[] rawImage, int imageW, int imageH, bool bTop, bool bBtm, double fCutRatio, byte value)
        {
            // position exception
            if (fCutRatio <= 0 || fCutRatio > 1.0) return;

            int nCutPos = 0;

            if (bTop)
            {
                nCutPos = Convert.ToInt32(imageH * fCutRatio);

                Parallel.For(0, imageW, x =>
                {
                    for (int y = 0; y < nCutPos; y++)
                    {
                        rawImage[y * imageW + x] = value;
                    }
                });
            }

            if (bBtm)
            {
                nCutPos = Convert.ToInt32(imageH * (1 - fCutRatio));

                Parallel.For(0, imageW, x =>
                {
                    for (int y = nCutPos; y < imageH; y++)
                    {
                        rawImage[y * imageW + x] = value;
                    }
                });
            }
        }
        public static void HC_FILL_CircleRegion(byte[] rawImage, int imageW, int imageH, bool bMajorAxisBased, byte value)
        {
            // radius = major langth / 2.0;
            // radius = minor length / 2.0;
            int nRadius = 0;

            if (bMajorAxisBased == true)
                nRadius = (int)((Math.Max(imageW, imageH) / 2.0));
            else
                nRadius = (int)((Math.Min(imageW, imageH) / 2.0));

            int nCenterX = imageW / 2;
            int nCenterY = imageH / 2;

            // get center
            PointF ptCenter = new PointF(nCenterX, nCenterX);

            //for (int y = 0; y < imageH; y++)
            Parallel.For(0, imageH, y =>
            {
                for (int x = 0; x < imageW; x++)
                {
                    CLine line = new CLine(x, y, nCenterX, nCenterY);
                    if (line.LENGTH > nRadius)
                    {
                        rawImage[y * imageW + x] = value;
                    }
                }
            });
        }

        public static void HC_FILL_RectRegion(byte[] rawImage, int imageW, int imageH, RectangleF rc, byte value)
        {
            int limitX = Convert.ToInt32(CRect.GetRB(rc).X);
            int limitY = Convert.ToInt32(CRect.GetRB(rc).Y);

            Parallel.For((int)rc.Y, limitY, y =>
            {
                for (int x = (int)rc.X; x < limitX; x++)
                {
                    rawImage[y * imageW + x] = value;
                }
            });
        }

        /// <summary>
        /// 비율에 따른 영상 커팅 
        /// </summary>
        /// <param name="rawImage"></param>
        /// <param name="imageW"></param>
        /// <param name="imageH"></param>
        /// <param name="fFactor"></param>
        /// <returns></returns>
        public static byte[] HC_FILL_WeightFiller(byte[] rawImage, int imageW, int imageH, double fFactor)
        {
            Array.Copy(rawImage, rawImage, rawImage.Length);

            int nMax = rawImage.Max();

            double fWeight = fFactor / 100.0;
            int threshold = Convert.ToInt32(nMax - (nMax * fWeight));

            Parallel.For(0, rawImage.Length, i =>
            {
                if (rawImage[i] < threshold)
                {
                    rawImage[i] = 0;
                }
            });
            return rawImage;
        }
    }
}
