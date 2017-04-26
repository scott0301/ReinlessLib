using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Drawing;
using System.Drawing.Imaging;

namespace ReinlessLib
{
    public static partial class CReinlessLib
    {

        #region IMAGE SAVE
        public static byte[] /***/LoadImage(string strPath, ref int imageW, ref int imageH)
        {
            Bitmap bmp = (Bitmap)Image.FromFile(strPath);

            byte[] rawImage = HC_CONV_Bmp2Byte(bmp, ref imageW, ref imageH);
            return rawImage;

        }
        public static void /*****/SaveImage(byte[] rawImage, int imageW, int imageH, string strPath)
        {
            Bitmap bmp = (Bitmap)HC_CONV_Byte2Bmp(rawImage, imageW, imageH);

            bmp.Save(strPath);
        }
        public static void /*****/SaveImage(double[] fImage, int imageW, int imageH, string strPath)
        {
            byte[] rawImage = HC_CONV_Double2Byte(fImage);

            Bitmap bmp = (Bitmap)HC_CONV_Byte2Bmp(rawImage, imageW, imageH);
            bmp.Save(strPath);
        }
        #endregion

        

        public static void HC_ARRAY_Dump(string strPath, double[] fArray, int w, int h)
        {
            string strBody = string.Empty;

            if (fArray.Length < w * h) return;

            for (int y = 0; y < h; y++)
            {
                for (int x = 0; x < w; x++)
                {
                    strBody += fArray[y * w + x].ToString() + ",";
                }
                strBody += System.Environment.NewLine;
            }

            try
            {

                System.IO.File.WriteAllText(strPath, strBody);
            }
            catch { }
        }
        public static void HC_ARRAY_Dump(string strPath, int[] nArray)
        {

            string strBody = string.Empty;

            for (int i = 0; i < nArray.Length; i++)
            {
                strBody += nArray[i].ToString() + ",";
            }

            try
            {
                System.IO.File.WriteAllText(strPath, strBody);
            }
            catch { }
        }
        public static void HC_ARRAY_Dump(string strPath, byte[] byteArray, int w, int h)
        {
            string strBody = string.Empty;

            for (int y = 0; y < h; y++)
            {
                for (int x = 0; x < w; x++)
                {
                    strBody += byteArray[y * w + x].ToString() + ",";
                }
                strBody += System.Environment.NewLine;
            }
            System.IO.File.WriteAllText(strPath, strBody);
        }

        /// <summary> 170426
        /// Set the element values : all clear  = object type 
        /// </summary>
        /// <param name="array"></param>
        /// <param name="fvalue"></param>
        public static void HC_ARRAY_SetValue(ref object array, double fvalue)
        {
            if (array.GetType() == new double[0].GetType())
            {
                double[] fArray = (double[])array;
                array = Enumerable.Repeat(fvalue, fArray.Length);
            }
            else if (array.GetType() == new float[0].GetType())
            {
                float[] fArray = (float[])array;
                array = Enumerable.Repeat((float)fvalue, fArray.Length);
            }
            else if (array.GetType() == new int[0].GetType())
            {
                int[] nArray = (int[])array;
                array = Enumerable.Repeat((int)fvalue, nArray.Length);
            }
            else if (array.GetType() == new byte[0].GetType())
            {
                byte[] byteArray = (byte[])array;
                byte c = fvalue < 0 ? (byte)0 : fvalue > 255 ? (byte)255 : (byte)fvalue;
                array = Enumerable.Repeat(c, byteArray.Length);
            }

        }
        
        /// <summary> 170426
        /// Get the Max element Position From  the Input array = object type 
        /// </summary>
        /// <param name="array"></param>
        /// <returns></returns>
        public static int HC_ARRAY_GetMaxElementPosition(object array)
        {
            int nIndex = 0;

            if (array.GetType() == new double[0].GetType())
            {
                double[] arrConv = (double[])array;
                double max = arrConv.Max();
                nIndex = Array.IndexOf(arrConv, max);
            }
            else if (array.GetType() == new float[0].GetType())
            {
                float[] arrConv = (float[])array;
                float max = arrConv.Max();
                nIndex = Array.IndexOf(arrConv, max);
            }
            else if (array.GetType() == new int[0].GetType())
            {
                int[] arrConv = (int[])array;
                int max = arrConv.Max();
                nIndex = Array.IndexOf(arrConv, max);
            }
            else if (array.GetType() == new byte[0].GetType())
            {
                byte[] arrConv = (byte[])array;
                byte max = arrConv.Max();
                nIndex = Array.IndexOf(arrConv, max);
            }
            return nIndex;
        }
        /// <summary> 170426
        /// Get the Max element Position From  the Input array = object type 
        /// </summary>
        /// <param name="array"></param>
        /// <returns></returns>
        public static int HC_ARRAY_GetMinElementPosition(object array)
        {
            int nIndex = 0;

            if (array.GetType() == new double[0].GetType())
            {
                double[] arrConv = (double[])array;
                double min = arrConv.Min();
                nIndex = Array.IndexOf(arrConv, min);
            }
            else if (array.GetType() == new float[0].GetType())
            {
                float[] arrConv = (float[])array;
                float min = arrConv.Min();
                nIndex = Array.IndexOf(arrConv, min);
            }
            else if (array.GetType() == new int[0].GetType())
            {
                int[] arrConv = (int[])array;
                int min = arrConv.Min();
                nIndex = Array.IndexOf(arrConv, min);
            }
            else if (array.GetType() == new byte[0].GetType())
            {
                byte[] arrConv = (byte[])array;
                byte min = arrConv.Min();
                nIndex = Array.IndexOf(arrConv, min);
            }
            return nIndex;
        }

        public static int HC_ARRAY_GetMatchedCount(object array, double fValue)
        {
            int nCount = 0;

            if (array.GetType() == new double[0].GetType())
            {
                double[] arrConv = (double[])array;
                nCount = Array.FindAll(arrConv, element => element == fValue).Length;
            }
            else if (array.GetType() == new float[0].GetType())
            {
                float[] arrConv = (float[])array;
                nCount = Array.FindAll(arrConv, element => element == (float)fValue).Length;
            }
            else if (array.GetType() == new int[0].GetType())
            {
                int[] arrConv = (int[])array;
                nCount = Array.FindAll(arrConv, element => element == (int)fValue).Length;
            }
            else if (array.GetType() == new byte[0].GetType())
            {
                byte[] arrConv = (byte[])array;
            }
            return nCount;
        }

        static int SortDoubleInt(KeyValuePair<double, int> a, KeyValuePair<double, int> b)
        {
            return b.Key.CompareTo(a.Key);
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
    }
}
