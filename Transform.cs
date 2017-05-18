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
        //******************************************************************************************
        // Condition
        //******************************************************************************************

        public static byte[] HC_TRANS_Reverse(byte[] rawImage, int imageW, int imageH)
        {
            byte[] newRaw = new byte[imageW * imageH];
            Parallel.For(0, imageW, x => { for (int y = 0; y < imageH; y++) { newRaw[y * imageW + x] = (byte)(255 - rawImage[y * imageW + x]); } });
            return newRaw;
        }
        public static byte[] HC_TRANS_Contrast(byte[] rawImage, int imageW, int imageH, int nValue)
        {
            double fContrastLevel = Math.Pow((100 + nValue) / 100.0, 2.0);

            byte[] newRes = new byte[imageW * imageH];

            for (int i = 0; i < newRes.Length; i++)
            {
                double newValue = (((((rawImage[i] / 255.0) - 0.5) * fContrastLevel) + 0.5) * 255.0);

                newRes[i] = newValue > 255 ? (byte)255 : newValue > 0 ? (byte)newValue : (byte)0;

            }
            return newRes;
        }
        public static byte[] HC_TRANS_ContrastLUT(byte[] rawImage, int imageW, int imageH, int nValue)
        {
            double fSigma = 0;

            if (nValue <= 0)
            {
                fSigma = 1.0 + (nValue / 256.0);
            }
            else
            {
                fSigma = 256.0 / Math.Pow(2, Math.Log(257 - nValue, 2));
            }

            byte[] LUT = new byte[256];

            for (int i = 0; i < 256; i++)
            {
                if ((fSigma * (i - 127) + 127) > 255)
                {
                    LUT[i] = 255;
                }
                else if ((fSigma * (i - 127) + 127) < 0)
                {
                    LUT[i] = 0;
                }
                else
                {
                    double f = fSigma * ((double)i - 127.0) + 127.0;

                    LUT[i] = EnsureByte(f);
                }
            }


            for (int i = 0; i < rawImage.Length; i++)
            {
                rawImage[i] = LUT[rawImage[i]];
            }

            return rawImage;
        }

        // 170428 sigmoidar contarst correction
        public static byte[] HC_TRANS_SIGMOID_Contrast(byte[] rawImage, int imageW, int imageH, double fCutoff, double fGain)
        {
            double[] fImage = CReinlessLib.HC_CONV_Byte2Double(rawImage);

            double min = fImage.Min();
            Parallel.For(0, fImage.Length, i => { fImage[i] -= min; });

            // devision max
            double max = fImage.Max();

            Parallel.For(0, fImage.Length, i => { fImage[i] /= max; });

            // case 2 : gamma 
            Parallel.For(0, fImage.Length, i =>
            {
                double fValue = (fCutoff - fImage[i]) * fGain;
                fImage[i] = 1.0 / (1.0 + Math.Exp(fValue));
            });

            return CReinlessLib.HC_CONV_GetNormalizedImage(fImage);

        }
        // 170428 sigmoidal gamma correction
        public static byte[] HC_TRANS_SIGMOID_Gamma(byte[] rawImage, int imageW, int imageH, double fCutoff, double fGain)
        {
            double[] fImage = CReinlessLib.HC_CONV_Byte2Double(rawImage);

            double min = fImage.Min();
            Parallel.For(0, fImage.Length, i => { fImage[i] -= min; });

            // devision max
            double max = fImage.Max();

            Parallel.For(0, fImage.Length, i => { fImage[i] /= max; });

            // case 2 : gamma 
            Parallel.For(0, fImage.Length, i =>
            {
                double fValue = (fCutoff - fImage[i]) * fGain;
                fImage[i] = Math.Pow(fImage[i], 1.0 / fGain);
            });

            return CReinlessLib.HC_CONV_GetNormalizedImage(fImage);

        }

        public static byte[] HC_TRANS_Brightness(byte[] rawImage, int imageW, int imageH, int nIncrement)
        {
            byte[] newRaw = new byte[imageW * imageH];

            Parallel.For(0, rawImage.Length, i =>
            {
                int nValue = rawImage[i] + nIncrement;
                nValue = nValue > 255 ? 255 : nValue < 0 ? 0 : nValue;
                newRaw[i] = (byte)nValue;
            });

            return newRaw;
        }
        public static byte[] HC_TRANS_QUANTIZATION(byte[] rawimage, int imageW, int imageH, int nLevel)
        {
            int nBit = 0;

            switch (nLevel)
            {
                case 128: nBit = 1; break;
                case 64: nBit = 2; break;
                case 32: nBit = 3; break;
                case 16: nBit = 4; break;
                case 8: nBit = 5; break;
                case 4: nBit = 6; break;
                case 2: nBit = 7; break;
                default: nBit = 0; break;
            }

            byte[] newRaw = new byte[imageW * imageH];
             

            for (int y = 0; y < imageH; y++)
            {
                for (int x = 0; x < imageW; x++)
                {
                    int nPixelValue = rawimage[y * imageW + x];
                    nPixelValue = nPixelValue >> nBit;
                    nPixelValue = nPixelValue << nBit;
                    newRaw[y * imageW + x] = (byte)nPixelValue;
                }
            }

            return newRaw;
        }

        public static byte[]/****/HC_TRANS_FlipX(byte[] rawImage, int imageW, int imageH)
        {
            byte[] flipImage = new byte[imageW * imageH];
            Parallel.For(0, imageH, y => { for (int x = 0, revX = imageW - 1; x < imageW; x++, revX--) { flipImage[y * imageW + revX] = rawImage[y * imageW + x]; } });
            return flipImage;
        }
        public static byte[]/****/HC_TRANS_FlipY(byte[] rawImage, int imageW, int imageH)
        {
            byte[] flipImage = new byte[imageW * imageH];
            Parallel.For(0, imageW, x => { for (int y = 0, revY = imageH - 1; y < imageH; y++, revY--) { flipImage[revY * imageW + x] = rawImage[y * imageW + x]; } });
            return flipImage;
        }

        // black average value ==> types of mosaic
        public static byte[] /****/HC_TRANS_Uniformity(byte[] rawImage, int imageW, int imageH)
        {
            byte[] newRaw = new byte[imageW * imageH];

            Parallel.For(0, 15, loop =>
            {
                for (int block = 1; block < 16; block++)
                {
                    for (int y = 0; y < imageH; y += block)
                    {
                        for (int x = 0; x < imageW; x += block)
                        {
                            double fx = 0;
                            double fcount = 0;

                            for (int yy = y; yy < y + block; yy++)
                            {
                                for (int xx = x; xx < x + block; xx++)
                                {
                                    if (yy >= imageH) continue;
                                    if (xx >= imageW) continue;

                                    fx += rawImage[yy * imageW + xx];
                                    fcount++;
                                }
                            }

                            fx /= fcount;

                            for (int yy = y; yy < y + block; yy++)
                            {
                                for (int xx = x; xx < x + block; xx++)
                                {
                                    if (yy >= imageH) continue;
                                    if (xx >= imageW) continue;

                                    newRaw[yy * imageW + xx] = Convert.ToByte(fx);
                                }
                            }

                        }
                    }
                }
            });

            return newRaw;
        }

        public static byte[] HC_TRANS_GradientImage(byte[] rawImage, int imageW, int imageH)
        {
            double[] fImage = new double[imageW * imageH];

            //for (int y = 1; y < imageH - 1; y++)
            Parallel.For(1, imageH - 1, y =>
            {
                for (int x = 1; x < imageW - 1; x++)
                {
                    double dx = (rawImage[y * imageW + x + 1] - rawImage[y * imageW + x - 1]) * (rawImage[y * imageW + x + 1] - rawImage[y * imageW + x - 1]);
                    double dy = (rawImage[(y + 1) * imageW + x] - rawImage[(y - 1) * imageW + x]) * (rawImage[(y + 1) * imageW + x] - rawImage[(y - 1) * imageW + x]);
                    double fValue = Math.Sqrt(dx + dy);//Math.Sqrt(dx + dy);

                    if (double.IsNaN(fValue)) fValue = 0;

                    fImage[y * imageW + x] = fValue;
                }
            });
            return HC_CONV_GetNormalizedImage(fImage);
        }
        public static byte[] HC_TRANS_GradientX(byte[] rawImage, int imageW, int imageH)
        {
            double[] fImage = new double[imageW * imageH];

            //for (int y = 1; y < imageH - 1; y++)
            Parallel.For(1, imageH - 1, y =>
            {
                for (int x = 1; x < imageW - 1; x++)
                {
                    double dx = (rawImage[y * imageW + x + 1] - rawImage[y * imageW + x - 1]) * (rawImage[y * imageW + x + 1] - rawImage[y * imageW + x - 1]);
                    double dy = (rawImage[(y + 1) * imageW + x] - rawImage[(y - 1) * imageW + x]) * (rawImage[(y + 1) * imageW + x] - rawImage[(y - 1) * imageW + x]);
                    double fValue = Math.Sqrt(dx * dx);//Math.Sqrt(dx + dy);

                    if (double.IsNaN(fValue)) fValue = 0;

                    fImage[y * imageW + x] = fValue;
                }
            });
            return HC_CONV_GetNormalizedImage(fImage);
        }
        public static byte[] HC_TRANS_GradientY(byte[] rawImage, int imageW, int imageH)
        {
            double[] fImage = new double[imageW * imageH];

            //for (int y = 1; y < imageH - 1; y++)
            Parallel.For(1, imageH - 1, y =>
            {
                for (int x = 1; x < imageW - 1; x++)
                {
                    double dx = (rawImage[y * imageW + x + 1] - rawImage[y * imageW + x - 1]) * (rawImage[y * imageW + x + 1] - rawImage[y * imageW + x - 1]);
                    double dy = (rawImage[(y + 1) * imageW + x] - rawImage[(y - 1) * imageW + x]) * (rawImage[(y + 1) * imageW + x] - rawImage[(y - 1) * imageW + x]);
                    double fValue = Math.Sqrt(dy + dy);//Math.Sqrt(dx + dy);

                    if (double.IsNaN(fValue)) fValue = 0;

                    fImage[y * imageW + x] = fValue;
                }
            });
            return HC_CONV_GetNormalizedImage(fImage);
        }

        //******************************************************************************************
        // Resize
        //******************************************************************************************

        #region RESIZE
        public static byte[] /***/HC_TRANS_RESIZE_NearestNeibor(byte[] rawImage, int imageW, int imageH, int resizedW, int resizedH)
        {
            byte[] scaledImage = new byte[resizedW * resizedH];

            // get sub sample image
            for (int y = 0; y < imageH; y++)
            {
                for (int x = 0; x < imageW; x++)
                {
                    scaledImage[(y * 2) * resizedW + (x * 2)] = rawImage[y * imageW + x];
                }
            }

            // horizontal filling
            for (int y = 0; y < resizedH - 1; y += 2)
            {
                for (int x = 1; x < resizedW - 1; x += 2)
                {
                    int val1 = scaledImage[y * resizedW + x - 1];
                    int val2 = scaledImage[y * resizedW + x + 1];
                    int res = Convert.ToInt32((val1 + val2) / 2.0);
                    if (res > 255) res = 255;
                    scaledImage[y * resizedW + x] = (byte)res;
                }
            }

            //vertical filling
            for (int y = 1; y < resizedH - 1; y += 2)
            {
                for (int x = 0; x < resizedW - 1; x += 2)
                {
                    int val1 = scaledImage[(y - 1) * resizedW + x];
                    int val2 = scaledImage[(y + 1) * resizedW + x];

                    int res = Convert.ToInt32((val1 + val2) / 2.0);
                    if (res > 255) res = 255;
                    scaledImage[y * resizedW + x] = (byte)res;
                }
            }

            //vertical filling
            for (int y = 1; y < resizedH - 1; y += 2)
            {
                for (int x = 1; x < resizedW - 1; x += 2)
                {
                    int val1 = scaledImage[(y + 0) * resizedW + x - 1];
                    int val2 = scaledImage[(y + 0) * resizedW + x + 1];
                    int val3 = scaledImage[(y - 1) * resizedW + x];
                    int val4 = scaledImage[(y + 1) * resizedW + x];

                    int res = Convert.ToInt32((val1 + val2 + val3 + val4) / 4.0);
                    if (res > 255) res = 255;
                    scaledImage[y * resizedW + x] = (byte)res;
                }
            }
            return scaledImage;
        }
        public static byte[] /***/HC_TRANS_RESIZE_Bicubic(byte[] rawImage, int imageW, int imageH, int resizedW, int resizedH)
        {
            byte[] temp = new byte[resizedW * resizedH];
            int index;

            float x_ratio = ((float)(imageW - 1)) / resizedW;
            float y_ratio = ((float)(imageH - 1)) / resizedH;

            float x_diff, y_diff;
            int offset = 0;
            for (int i = 0; i < resizedH; i++)
            //Parallel.For(0, resizedH, i =>
            {
                for (int j = 0; j < resizedW; j++)
                {
                    int x = (int)(x_ratio * j);
                    int y = (int)(y_ratio * i);

                    x_diff = (x_ratio * j) - x;
                    y_diff = (y_ratio * i) - y;
                    index = y * imageW + x;

                    // range is 0 to 255 thus bitwise AND with 0xff
                    int A = rawImage[index] & 0xff;
                    int B = rawImage[index + 1] & 0xff;
                    int C = rawImage[index + imageW] & 0xff;
                    int D = rawImage[index + imageW + 1] & 0xff;

                    // Y = A(1-w)(1-h) + B(w)(1-h) + C(h)(1-w) + Dwh
                    int gray = (int)(A * (1 - x_diff) * (1 - y_diff) + B * (x_diff) * (1 - y_diff) + C * (y_diff) * (1 - x_diff) + D * (x_diff * y_diff));

                    temp[offset++] = (byte)gray;
                }
            }//);
            return temp;
        }
        #endregion

        #region CROP
        public static byte[] HC_CropImage(byte[] rawInput, int imageW, int imageH, Rectangle rc)
        {
            int x = rc.X;
            int y = rc.Y;
            int cropW = rc.Width;
            int cropH = rc.Height;

            return HC_CropImage(rawInput, imageW, imageH, x, y, cropW, cropH);

        }
        public static byte[] HC_CropImage(byte[] rawInput, int imageW, int imageH, int ptX, int ptY, int cropW, int cropH)
        {
            byte[] rawCrop = new byte[cropW * cropH];


            for (int y = ptY, copyLine = 0; y < ptY + cropH; y++)
            {
                Buffer.BlockCopy(rawInput, y * imageW + ptX, rawCrop, cropW * copyLine++, cropW);
            }

            return rawCrop;

        }
        public static byte[] HC_CropImage(byte[] rawInput, int imageW, int imageH, RectangleF rc)
        {
            int nLength = Convert.ToInt32(rc.Width * rc.Height);
            int toHeight = Convert.ToInt32(rc.Y + rc.Height);
            int toWidth = Convert.ToInt32(rc.Width);
            int px = Convert.ToInt32(rc.X);


            byte[] rawCrop = new byte[nLength];

            for (int y = (int)rc.Y, copyLine = 0; y < toHeight; y++)
            {
                Buffer.BlockCopy(rawInput, y * imageW + px, rawCrop, toWidth * copyLine++, toWidth);
            }

            return rawCrop;
        }

        // Polar croodinate cropping
        // 입력된 rectangle의 polar croodinate에 대한 cropping 결과를 return --> pixelation event occured. --> interpolation 적용
        public static byte[] HC_CropImage_Polar(byte[] rawInput, int imageW, int imageH, RectangleF rc)
        {
            int nCX = Convert.ToInt32(rc.Width / 2.0);
            int nCY = Convert.ToInt32(rc.Height / 2.0);
            int nRadius = Math.Max(nCX, nCY);

            nCX += (int)rc.X;
            nCY += (int)rc.Y;

            byte[] rawPolar = new byte[360 * nRadius];

            for (int na = 0; na < 360; na++)
            {
                int nPolarY = nRadius;
                for (int nr = 0; nr < nRadius; nr++)
                {
                    double fDegree = ((na - 90.0) * Math.PI / 180.0);

                    double x = nCX + (nr * Math.Cos(fDegree));
                    double y = nCY + (nr * Math.Sin(fDegree));

                    if (x < 0 || y < 0 || x >= imageW || y >= imageH)
                    {
                        continue;
                    }
                    rawPolar[--nPolarY * 360 + na] = rawInput[(int)y * imageW + (int)x];
                }
            }
            return rawPolar;
        }
        // interpolated polar transfrom --> plxelation occured --> enhanced version 170518 updated
        public static byte[] HC_CropImage_Interpolated_Polar(byte[] rawInput, int imageW, int imageH, RectangleF rc)
        {
            int nCX = Convert.ToInt32(rc.Width / 2.0);
            int nCY = Convert.ToInt32(rc.Height / 2.0);
            int nRadius = Math.Max(nCX, nCY);

            nCX += (int)rc.X;
            nCY += (int)rc.Y;

            byte[] rawPolar = new byte[360 * nRadius];

            for (int na = 0; na < 360; na++)
            {
                int nPolarY = nRadius;
                for (int nr = 0; nr < nRadius; nr++)
                {
                    double fDegree = ((na - 90.0) * Math.PI / 180.0);

                    double x = nCX + ((double)nr * Math.Cos(fDegree));
                    double y = nCY + ((double)nr * Math.Sin(fDegree));

                    if (x < 0 || y < 0 || x >= imageW || y >= imageH) { continue; }

                    int x1 = (int)Math.Floor(x);
                    int x2 = (int)Math.Ceiling(x);
                    int y1 = (int)Math.Floor(y);
                    int y2 = (int)Math.Ceiling(y);

                    int q11 = rawInput[y1 * imageW + x1];
                    int q12 = rawInput[y2 * imageW + x1];
                    int q21 = rawInput[y1 * imageW + x2];
                    int q22 = rawInput[y2 * imageW + x2];

                    double fInterplated = GetInterPolatedValue(x, y, x1, x2, y1, y2, q11, q12, q21, q22);

                    byte valueOrg = rawInput[(int)y * imageW + (int)x];

                    byte value = fInterplated < 0 ? (byte)0 : fInterplated > 255 ? (byte)255 : (byte)fInterplated;

                    if (double.IsNaN(fInterplated) == true) value = valueOrg;

                    rawPolar[--nPolarY * 360 + na] = value;
                }
            }
            return rawPolar;
        }


        // rotated rectangle Cropping
        // Rotated Recangle의 Unroated Rectangle을 Cropping 하여 반환
        public static byte[] HC_CropImage_Rotate(byte[] rawInput, int imageW, int imageH, RectangleF rc, PointF ptGravity, float fAngle)
        {
            List<byte> listRot = new List<Byte>();

            int fromX = (int)(rc.X);
            int fromY = (int)(rc.Y);

            for (int y = fromY; y < (int)fromY + (int)rc.Height; y++)
            {
                for (int x = fromX; x < (int)fromX + (int)rc.Width; x++)
                {
                    PointF ptRot = _RotatePointByGravity(new PointF(x, y), ptGravity, fAngle);


                    double cx = ptRot.X;
                    double cy = ptRot.Y;
                    int x1 = (int)Math.Floor(cx);
                    int x2 = (int)Math.Ceiling(cx);
                    int y1 = (int)Math.Floor(cy);
                    int y2 = (int)Math.Ceiling(cy);

                    int q11 = rawInput[y1 * imageW + x1];
                    int q12 = rawInput[y2 * imageW + x1];
                    int q21 = rawInput[y1 * imageW + x2];
                    int q22 = rawInput[y2 * imageW + x2];

                    double fInterplated = GetInterPolatedValue(cx, cy, x1, x2, y1, y2, q11, q12, q21, q22);

                    byte value = fInterplated < 0 ? (byte)0 : fInterplated > 255 ? (byte)255 : (byte)fInterplated;
                    listRot.Add(value);
                }
            }

            byte[] rawRes = new byte[(int)rc.Width * (int)rc.Height];
            byte[] rawOut = listRot.ToArray();
            return rawOut;
        }
        public static byte[] HC_CropImage_Overlap(byte[] rawInput, int imageW, int imageH, byte[] cropImage, int cropW, int cropH, RectangleF rc)
        {
            int posX = Convert.ToInt32(rc.X);
            int posY = Convert.ToInt32(rc.Y);
            Parallel.For(0, cropH, y => { Buffer.BlockCopy(cropImage, y * cropW, rawInput, (posY + y) * imageW + posX, cropW); });
            return rawInput;
        }
        #endregion  

        #region Rotaiton
        public static Point _GetRotatePos(double X, double Y, double fAngle, int nOffsetX, int nOffsetY)
        {
            int nRotatedX = (int)(Math.Round((X - nOffsetX) * Math.Cos(fAngle) - (Y - nOffsetY) * Math.Sin(fAngle))) + nOffsetX;
            int nRotatedY = (int)(Math.Round((X - nOffsetY) * Math.Sin(fAngle) + (Y - nOffsetY) * Math.Cos(fAngle))) + nOffsetY;

            return new Point(nRotatedX, nRotatedY);
        }
        private static PointF _RotatePointByGravity(PointF ptTarget, PointF ptGravity, double fAngle)
        {
            //x' = (x-a) * cosR - (y-b)sinR + a
            //y' = (x-a) * sinR + (y-b)cosR + b
            fAngle = fAngle * Math.PI / 180.0;

            PointF ptRotated = new PointF(0, 0);

            ptRotated.X = (float)(((ptTarget.X - ptGravity.X) * Math.Cos(fAngle) - (ptTarget.Y - ptGravity.Y) * Math.Sin(fAngle)) + ptGravity.X);
            ptRotated.Y = (float)(((ptTarget.X - ptGravity.X) * Math.Sin(fAngle) + (ptTarget.Y - ptGravity.Y) * Math.Cos(fAngle)) + ptGravity.Y);

            return ptRotated;
        }
        public static double GetInterPolatedValue(double cx, double cy, double x1, double x2, double y1, double y2, double q11, double q12, double q21, double q22)
        {
            double r1 = (((x2 - cx) / (x2 - x1)) * q11) + (((cx - x1) / (x2 - x1)) * q21);
            double r2 = (((x2 - cx) / (x2 - x1)) * q12) + (((cx - x1) / (x2 - x1)) * q22);
            double pvalue = (((y2 - cy) / (y2 - y1)) * r1) + (((cy - y1) / (y2 - y1)) * r2);
            return pvalue;
        }
        
        public static byte[] HC_TRANS_RotationBasic(byte[] rawImage, int imageW, int imageH, double fAngle)
        {
            byte[] rotatedImage = new byte[imageW * imageH];

            double fDegree = fAngle * Math.PI / 180.0;

            int nOffsetX = imageW / 2;
            int nOffsetY = imageH / 2;

            for (int y = 0; y < imageH; y++)
            {
                for (int x = 0; x < imageW; x++)
                {
                    Point ptRotated = _GetRotatePos(x, y, fDegree, nOffsetX, nOffsetY);

                    if (ptRotated.X > 0 && ptRotated.X < imageW && ptRotated.Y > 0 && ptRotated.Y < imageH)
                    {
                        byte pixelValue = rawImage[y * imageW + x];
                        rotatedImage[ptRotated.Y * imageW + ptRotated.X] = pixelValue;
                    }

                }
            }

            return rotatedImage;
        }
        public static byte[] HC_TRANS_RotationInterpolated(byte[] rawImage, int imageW, int imageH, double fAngle)
        {
            byte[] rotatedImage = new byte[imageW * imageH];

            // rotating clockwise, so it's negative relative to Cartesian quadrants
            double fDegree = fAngle * Math.PI / 180.0;

            int nCX = imageW / 2;
            int nCY = imageH / 2;

            // assigning pixels of destination image from source image
            // with bilinear interpolation
            //for (int y = 0; y < imageH; y++)
            Parallel.For(0, imageH, y =>
            {
                int vy = nCY - y;

                for (int x = 0; x < imageW; x++)
                //Parallel.For(0, imageW, x =>
                {
                    // convert raster to Cartesian
                    int vx = x - nCX;

                    double fPolarAngle = 0.0;
                    if (vx == 0 && vy == 0)
                    {
                        // this is center case, no need to rotate
                        rotatedImage[y * imageW + x] = rawImage[y * imageW + x];
                    }
                    else // if this is not the center point
                    {
                        if (vx == 0)
                        {
                            if (vy < 0) { fPolarAngle = 1.5 * Math.PI; }
                            else { fPolarAngle = 0.5 * Math.PI; }
                        }
                        else
                        {
                            fPolarAngle = Math.Atan2((double)vy, (double)vx);
                        }

                        // the crucial rotation part
                        // "reverse" rotate, so minus instead of plus
                        fPolarAngle -= fDegree;

                        // convert polar to Cartesian
                        double fDistance = Math.Sqrt(vx * vx + vy * vy);
                        double fTrueX = fDistance * Math.Cos(fPolarAngle);
                        double fTrueY = fDistance * Math.Sin(fPolarAngle);

                        // convert Cartesian to raster
                        fTrueX = fTrueX + (double)nCX;
                        fTrueY = (double)nCY - fTrueY;

                        int iFloorX = (int)(Math.Floor(fTrueX));
                        int iFloorY = (int)(Math.Floor(fTrueY));
                        int iCeilingX = (int)(Math.Ceiling(fTrueX));
                        int iCeilingY = (int)(Math.Ceiling(fTrueY));

                        // check bounds
                        if (iFloorX < 0 || iCeilingX < 0 || iFloorX >= imageW || iCeilingX >= imageW ||
                            iFloorY < 0 || iCeilingY < 0 || iFloorY >= imageH || iCeilingY >= imageH)
                        {
                            // this is boundary exception
                        }
                        else
                        {
                            double fDeltaX = fTrueX - (double)iFloorX;
                            double fDeltaY = fTrueY - (double)iFloorY;

                            int pixel_TL = rawImage[iFloorY * imageW + iFloorX]; ;
                            int pixel_TR = rawImage[iFloorY * imageW + iCeilingX];
                            int pixel_BL = rawImage[iCeilingY * imageW + iFloorX];
                            int pixel_BR = rawImage[iCeilingY * imageW + iCeilingX];

                            // linearly interpolate horizontally between top neighbours
                            // linearly interpolate horizontally between bottom neighbours
                            double fInterpolatedValueTop = (1 - fDeltaX) * pixel_TL + fDeltaX * pixel_TR;
                            double fInterpolatedBtm = (1 - fDeltaX) * pixel_BL + fDeltaX * pixel_BR;

                            // linearly interpolate vertically between top and bottom interpolated results
                            int nPixelValue = (int)(Math.Round((1 - fDeltaY) * fInterpolatedValueTop + fDeltaY * fInterpolatedBtm));

                            nPixelValue = nPixelValue < 0 ? 0 : nPixelValue;
                            nPixelValue = nPixelValue > 255 ? 255 : nPixelValue;

                            rotatedImage[y * imageW + x] = (byte)nPixelValue;
                           
                        }
                    }// else

                } // loop x
            }); // loop y

            return rotatedImage;
        }
        #endregion  

        public static byte[] /****/HC_TRANS_Cartesian2Polar_Raw(byte[] rawImage, int imageW, int imageH)
        {
            int nCenterX = imageW / 2;
            int nCenterY = imageH / 2;

            //int nRadius = Convert.ToInt32(Math.Sqrt((imageW * imageW) + (imageH * imageH)) / 2.0);
            int nRadius = Math.Min(nCenterX, nCenterY);

            byte[] newRaw = new byte[nRadius * 360];

            for (int nAngle = 0; nAngle < 360; nAngle++)
            {
                int polarY = nRadius;
                // Inner Contour of deformed circle can not detect when radius is assumed square.
                // +50 additive gap allow to detect deformed circle
                for (int moveRadisu = 0; moveRadisu < nRadius; moveRadisu++)
                {
                    double fDegree =  ((nAngle - 90) * 3.14 / 180);

                    int x = nCenterX + (int)(moveRadisu * Math.Cos(fDegree)); 
                    int y = nCenterY + (int)(moveRadisu * Math.Sin(fDegree));

                    if (x < 0 || y < 0 || x >= imageW || y >= imageH)
                    {
                        continue;
                    }

                    newRaw[(--polarY) * 360 + nAngle] = rawImage[y * imageW + x];

                } // radius loop
            } // angle loop

            return newRaw;
        }
        public static int /*******/HC_TRANS_GetPolarLength(Byte[] byteArray, int nWidth, int nHeight)
        {
            for (int x = nWidth - 1; x >= 0; x--)
            {
                int nSum = 0;
                for (int y = 0; y < nHeight; y++)
                {
                    nSum += byteArray[y * nWidth + x];
                }
                if (nSum == 0) return x;
            }
            return 0;
        }


        public static Point ShearXY(int x, int y, double shearX,double shearY,int offsetX,int offsetY)
        {
            Point result = new Point();

            result.X = (int)(Math.Round(x + shearX * y));
            result.X -= offsetX;

            result.Y = (int)(Math.Round(y + shearY * x));
            result.Y -= offsetY;

            return result;
        }
        public static byte[] HC_TRANS_Sheer(byte[] rawImage, int imageW, int imageH, double sX, double sY)
        {

            byte [] resImage = new byte[imageW*imageH];

            int offsetX = (int)Math.Round(imageW * sX / 2.0);
            int offsetY = (int)Math.Round(imageH * sY / 2.0);

            Rectangle imageBounds = new Rectangle(0, 0, imageW, imageH);
            Point sourcePoint = new Point();
            Point resultPoint = new Point();

            for (int y = 0; y < imageH; y++)
            {
                for (int x = 0; x < imageW; x++)
                {
                    int sourceXY = y * imageW + x;

                    sourcePoint.X = x;
                    sourcePoint.Y = y;

                    if (sourceXY >= 0 && sourceXY+ 1 < resImage.Length)
                    {
                        resultPoint = ShearXY(x, y, sX, sY, offsetX, offsetY);
                        int resultXY = resultPoint.X * imageW + resultPoint.Y;

                        if (imageBounds.Contains(resultPoint) && resultXY >= 0 && resultXY < resImage.Length)
                        {
                            
                            if( resultXY + 1 < resImage.Length)
                            {
                                resImage[resultXY + 1] = rawImage[sourceXY];
                            }

                            if (resultXY - 1 >= 0)
                            {
                                resImage[resultXY - 1] = rawImage[sourceXY];
                            }

                            if (resultXY  < resImage.Length)
                            {
                                resImage[resultXY] = rawImage[sourceXY];
                            }
                        }
                    }

                }
            }
            return resImage;

        }
        #region registration
        public static byte[] HC_TRANS_STITCHING_Hor(byte[] image1, int w1, int h1, byte[] image2, int w2, int h2, int nTransX, int nTransY, ref int mergedW, ref int mergedH)
        {
            mergedW = w1 + w2 + nTransX;
            mergedH = Math.Max(h1, h2) + Math.Abs(nTransY);

            byte[] imageMerged = new byte[mergedW * mergedH];

            Rectangle rcCrop1 = new Rectangle(0, 0, w1, h1);
            if (nTransY < 0) rcCrop1.Y = -nTransY; // right image - attach --> left image movedown

            imageMerged = CReinlessLib.HC_CropImage_Overlap(imageMerged, mergedW, mergedH, image1, w1, h1, rcCrop1);

            Rectangle rcCrop2 = new Rectangle(w1, 0, w2, h2);
            rcCrop2.X += nTransX;
            if (nTransY > 0) rcCrop2.Y = nTransY; // + attach --> right image y set 0

            imageMerged = CReinlessLib.HC_CropImage_Overlap(imageMerged, mergedW, mergedH, image2, w2, h2, rcCrop2);

            return imageMerged;
        }
        public static byte[] HC_TRANS_STITCHING_Ver(byte[] image1, int w1, int h1, byte[] image2, int w2, int h2, int nTransX, int nTransY, ref int mergedW, ref int mergedH)
        {
            mergedW = Math.Max(w1, w2) + Math.Abs(nTransX);
            mergedH = h1 + h2 + nTransY;

            byte[] imageMerged = new byte[mergedW * mergedH];

            Rectangle rcCrop1 = new Rectangle(0, 0, w1, h1);
            if (nTransX < 0) rcCrop1.X = -nTransX; // right image - attach --> left image moveright

            imageMerged = CReinlessLib.HC_CropImage_Overlap(imageMerged, mergedW, mergedH, image1, w1, h1, rcCrop1);

            Rectangle rcCrop2 = new Rectangle(0, h1, w2, h2);
            rcCrop2.Y += nTransY;
            /***/if (nTransX > 0) rcCrop2.X = nTransX; // + attach --> right image y set 0

            imageMerged = CReinlessLib.HC_CropImage_Overlap(imageMerged, mergedW, mergedH, image2, w2, h2, rcCrop2);

            return imageMerged;
        }
        public static byte[] HC_TRANS_MergingHor(byte[] image1, int w1, int h1, byte[] image2, int w2, int h2, int nGap)
        {
            DRect rcCrop1 = new DRect(w1 - nGap, 0, nGap, h1);
            DRect rcCrop2 = new DRect(0, 0, nGap, h2);

            byte[] cropImage1 = CReinlessLib.HC_CropImage(image1, w1, h1, rcCrop1.ToRectangle());
            byte[] cropImage2 = CReinlessLib.HC_CropImage(image2, w2, h2, rcCrop2.ToRectangle());

            int cropW = nGap;
            int cropH = h1;

            byte[] mergedImage = CReinlessLib.HC_TRANS_STITCHING_Hor(cropImage1, cropW, cropH, cropImage2, cropW, cropH, 0, 0, ref cropW, ref cropH);

            return mergedImage;
        }
        public static byte[] HC_TRANS_MergingVer(byte[] image1, int w1, int h1, byte[] image2, int w2, int h2, int nGap)
        {
            DRect rcCrop1 = new DRect(0, h1 - nGap, w1, nGap);
            DRect rcCrop2 = new DRect(0, 0, w2, nGap);

            byte[] cropImage1 = CReinlessLib.HC_CropImage(image1, w1, h1, rcCrop1.ToRectangle());
            byte[] cropImage2 = CReinlessLib.HC_CropImage(image2, w2, h2, rcCrop2.ToRectangle());

            int cropW = w1;
            int cropH = nGap;

            byte[] mergedImage = CReinlessLib.HC_TRANS_STITCHING_Ver(cropImage1, cropW, cropH, cropImage2, cropW, cropH, 0, 0, ref cropW, ref cropH);

            return mergedImage;
        }
        #endregion

        public static byte[] HC_TRANS_Stretch(byte[] rawImage, int imageW, int imageH)
        {
            int pixels = imageW * imageH;
            double constant = 255.0 / pixels;

            int[] arrCDF = HC_HISTO_GetHistogram(rawImage, imageW, imageH, 255);
            double[] fImage = new double[imageW * imageH];

            //Convert arrays to cumulative distribution frequency data
            for (int i = 1; i <= 255; i++)
            {
                arrCDF[i] = arrCDF[i] + arrCDF[i - 1];
            }

            for (int y = 0; y < imageH; y++)
            {
                for (int x = 0; x < imageW; x++)
                {
                    byte pixel =  rawImage[y * imageW + x];
                    double value = arrCDF[pixel] * constant;

                    fImage[y * imageW + x] = value;
                }
            }
            rawImage = HC_CONV_GetNormalizedImage(fImage);
            return rawImage;
        }
        public static byte[] HC_TRANS_Hough(byte[] rawImage, ref int imageW, ref int imageH)
        {
           double hough_h = ((Math.Sqrt(2.0) * (double)(imageH>imageW?imageH:imageW)) / 2.0);  

           int accH = Convert.ToInt32(hough_h * 2.0); // -r -> +r  
           int accW = Convert.ToInt32(180);

           double radian = Math.PI / 180;
            byte [] imageAcc = new byte[accW*accH];
    
            double cx = imageW/2;  
            double cy = imageH/2;

            for (int y = 0; y < imageH; y++)
            {
                for (int x = 0; x < imageW; x++)
                {
                    if (rawImage[y * imageW + x] > 250)
                    {
                        for (int t = 0; t < 180; t++)
                        {
                            // x * cos(t) + y * sin(t)
                            double r = ((x - cx) * Math.Cos(t * radian)) + ((y - cy) * Math.Sin(t * radian));

                            int nIndex = (int)((Math.Round(r + hough_h) * 180.0)) + t;
                            if (nIndex >= 0 && nIndex < imageAcc.Length)
                            {
                                imageAcc[nIndex]++;
                            }
                        }
                    }
                }
            }

            imageW = accW;
            imageH = accH;
            return imageAcc;
        }

        public static List<CLine> HC_TRANS_HoughGetLines(byte[] accImage, int accW, int accH, int imageW, int imageH, int nThreshold)
        {
            int cx = imageW / 2;
            int cy = imageH / 2;

            int acc_cy = accH / 2;

            double radian =  Math.PI / 180;

            List<CLine> list = new List<CLine>();

            for (int y = 0; y < accH; y++)
            {
                for (int x = 0; x < accW; x++)
                {
                    int maxima = accImage[y * accW + x];

                    if( maxima >= nThreshold ) 
                    {
                        //**********************************************
                        for (int yy = -4; yy <= 4; yy++)
                        {
                            for (int xx = -4; xx <= 4; xx++)
                            {
                                int nIndex = (y + yy) * accW + (x + xx);
                                if (nIndex > 0 && nIndex < accImage.Length)
                                {
                                    if ( accImage[nIndex] > maxima)
                                    {
                                        maxima = accImage[nIndex];
                                        xx = yy = 5;
                                    }
                                }

                            }// loop xx
                        } // loop yy
                        //**********************************************

                        if (maxima > accImage[y * accW + x]) continue;


                        if (x >= 45 && x <= 135)
                        {
                            //y = (r - x cos(t)) / sin(t)  
                            float x1 = 0;
                            float y1 = (float)(((y - acc_cy) - ((x1 - cx) * Math.Cos(x * radian))) / Math.Sin(x * radian) + cy);  
                            float x2 = imageW - 0;
                            float y2 = (float)(((y - acc_cy) - ((x2 - cx) * Math.Cos(x * radian))) / Math.Sin(x * radian) + cy);

                            CLine line = new CLine(x1, y1, x2, y2);
                            list.Add(line);
                        }
                        else
                        {
                            //x = (r - y sin(t)) / cos(t);  
                             float y1 = 0;
                             float x1 = (float)(((y - acc_cy) - ((y1 - cy) * Math.Sin(x * radian))) / Math.Cos(x * radian) + cx);  
                             float y2 = imageH - 0;
                             float x2 = (float)(((y - acc_cy) - ((y2 - cy) * Math.Sin(x * radian))) / Math.Cos(x * radian) + cx);
                             CLine line = new CLine(x1, y1, x2, y2);
                             list.Add(line);
                        }

                    }// if maxima
                } // loop x
            } // loop y

            return list;
        }

        public static byte[] HC_TRANS_DISTANCE(byte[] dataInput)
        {
            int arrLength = dataInput.Length;

            int k = 0;
            
            int[] v = new int[arrLength];
            double[] z = new double[arrLength + 1];
            double[] buf = new double[arrLength];

            v[0] = 0;
            z[0] = Double.NegativeInfinity;
            z[1] = Double.PositiveInfinity;

            double s = 0;

            for (int q = 1; q < arrLength; q++)
            {
                while (true)
                {
                    s = (((dataInput[q] + q * q) - (dataInput[v[k]] + v[k] * v[k])) / (2.0 * q - 2.0 * v[k]));

                    if (s <= z[k])
                    {
                        k--;
                    }
                    else
                    {
                        break;
                    }
                }

                k++;

                v[k] = q;
                z[k] = s;
                z[k + 1] = Double.PositiveInfinity;
            }

            k = 0;

            for (int q = 0; q < arrLength; q++)
            {
                while (z[k + 1] < q)
                {
                    k++;
                }

                buf[q] = ((q - v[k]) * (q - v[k]) + dataInput[v[k]]);
            }

            byte[] rawOut = new byte[arrLength];

            HC_CONV_Double2Byte(buf, rawOut);

            return rawOut;

        }// function 

        
        // 170119 
        public static List<CLine> HoughLines(byte [] rawImage, int imageW, int imageH, int nThresholHit, int nThreshold, double fThetaResolution, int nResultLimit )
        {
	        int nLengthDigonal = (int)(Math.Sqrt((double)(imageW*imageW + imageH*imageH)));  
	        int numRho = nLengthDigonal*2; // rho의 음수 영역을 위해 2를 곱함
	        int nThetaCount = Convert.ToInt32(180.0 / fThetaResolution);
	        int nTotalCount = numRho*nThetaCount; // rho와 theta 조합의 출현 횟수를 저장하는 공간

	        double []  LUT_Sin = new double[nThetaCount]; // sin 함수 룩업 테이블
	        double []  LUT_Cos = new double[nThetaCount]; // cos 함수 룩업 테이블

	        double toRad = Math.PI/nThetaCount;

            for (int nTheta = 0; nTheta < nThetaCount; nTheta++)
            {
                LUT_Sin[nTheta] = (double)Math.Sin(nTheta * toRad);
                LUT_Cos[nTheta] = (double)Math.Cos(nTheta * toRad);
            }

	        int [] arrTransformed = new int[nTotalCount];

            for (int y = 0; y < imageH; y++)
            {
                for (int x = 0; x < imageW; x++)
                {
                    if (rawImage[y * imageW + x] > nThreshold) // 경계선 픽셀
                    {
                        for (int nTheta = 0; nTheta < nThetaCount; nTheta++)
                        {
                            int rho = (int)(x * LUT_Sin[nTheta] + y * LUT_Cos[nTheta] + nLengthDigonal + 0.5);
                            arrTransformed[rho * nThetaCount + nTheta]++;
                        }
                    }
                }
            }

            List<int> listRho = new List<int>();
            List<double> listTheta = new List<double>();

	        // nThreshold을넘는 결과 저장
	        int nLine = 0;
            for (int i = 0; i < nTotalCount && nLine < nResultLimit; i++)
            {
                if (arrTransformed[i] > nThresholHit)
                {
                    int nIndex  = (int)(i / nThetaCount); // rho의 인덱스
                    listTheta.Add( (i - nIndex * nThetaCount) * fThetaResolution); //theta 의 인덱스
                    listRho.Add( nIndex - nLengthDigonal); // 음수 값이 차지하는 위치만큼 뺄셈
                    nLine++;
                }
            }

            List<CLine> list = new List<CLine>();

            for (int i = 0; i < listRho.Count; i++)
            {
                int nRho = listRho.ElementAt(i);
                double fTheta = listTheta.ElementAt(i);

                if (fTheta == 90) // 수직선
                {
                    CLine line = new CLine((float)nRho, 0, (float)nRho, imageH);
                    list.Add(line);
                }
                else
                {
                    int x1 = 0;
                    int y1 = (int)(nRho / Math.Cos(fTheta * Math.PI / 180.0) + 0.5);
                    int x2 = imageW;
                    int y2 = (int)((nRho - x2 * Math.Sin(fTheta * Math.PI / 180.0)) / Math.Cos(fTheta * Math.PI / 180.0) + 0.5);

                    CLine line = new CLine(x1, y1, x2, y2);
                    list.Add(line);
                }
            }

	        return list;
        }

        
    }
}
