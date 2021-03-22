using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

//
// 과제 1-2.
//
// 단위 변환된 파일들의 보간법(Interpolation)
// 파장범위(350~1000)
//
//2021.03.22 조계성

//3중 대각 행렬
namespace TridiagonalMatrix
{
    #region TRIDIAGONAL_MATRIX

    // Thomas algorithmn을 이용한 3중 대각 행렬 계산
    // https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    //
    public sealed class Tridiagonal
    {
        public double[] Solve(double[,] matrix, double[] d)
        {
            int rows = matrix.GetLength(0);
            int cols = matrix.GetLength(1);
            int len = d.Length;
            // 행, 열 개수가 동일한 정방 행렬일 때 동작
            if (rows == cols && rows == len)
            {

                double[] b = new double[rows];//0,0부터의 대각행렬
                double[] a = new double[rows];//1,0부터의 대각행렬
                double[] c = new double[rows];//0,1부터의 대각행렬
                //ax[i-1] + bx[i] + cx[i+1] = d 
                //[b1 c1        ][x1]   = [d1]
                //[a2 b2 c2     ][x2]   = [d2]
                //[   a3 b3 c3  ][x3]   = [d3]
                //[         cn-1][xn-1] = [dn-1]  
                //[        an bn][xn]   = [dn]

                // 3중 대각 행렬의 행과 열 요소 분리
                for (int i = 0; i < rows; i++)
                {
                    for (int j = 0; j < cols; j++)
                    {
                        if (i == j)             //b대각
                            b[i] = matrix[i, j];
                        else if (i == (j - 1))  //c대각
                            c[i] = matrix[i, j];
                        else if (i == (j + 1))  //a대각
                            a[i] = matrix[i, j];
                    }
                }
                try
                {
                    c[0] = c[0] / b[0];
                    d[0] = d[0] / b[0];

                    for (int i = 1; i < len - 1; i++)
                    {
                        c[i] = c[i] / (b[i] - a[i] * c[i - 1]);
                        d[i] = (d[i] - a[i] * d[i - 1]) / (b[i] - a[i] * c[i - 1]);
                    }
                    d[len - 1] = (d[len - 1] - a[len - 1] * d[len - 2]) / (b[len - 1] - a[len - 1] * c[len - 2]);

                    // back-substitution step
                    for (int i = (len - 1); i-- > 0;)
                    {
                        d[i] = d[i] - c[i] * d[i + 1];
                    }

                    return d;
                }
                catch (DivideByZeroException)   //b[0]가 0일때 발생하는 예외처리
                {
                    Console.WriteLine("b[0]=0");
                    return null;
                }
            }
            else
            {
                Console.WriteLine("입력한 행렬이 정방행렬이 아닙니다.");
                return null;
            }
        }
    }

    #endregion
}

namespace BasicInterpolation
{
    #region EXTENSIONS
    public static class Extensions
    {
        // returns a list sorted in ascending order 
        public static List<T> SortedList<T>(this T[] array)
        {
            List<T> l = array.ToList();
            l.Sort();
            return l;

        }

        // returns a difference between consecutive elements of an array
        public static double[] Diff(this double[] array)
        {
            int len = array.Length - 1;
            double[] diffsArray = new double[len];
            for (int i = 1; i <= len; i++)
            {
                diffsArray[i - 1] = array[i] - array[i - 1];
            }
            return diffsArray;
        }

        // scaled an array by another array of doubles
        public static double[] Scale(this double[] array, double[] scalor, bool mult = true)
        {
            int len = array.Length;
            double[] scaledArray = new double[len];

            if (mult)
            {
                for (int i = 0; i < len; i++)
                {
                    scaledArray[i] = array[i] * scalor[i];
                }
            }
            else
            {
                for (int i = 0; i < len; i++)
                {
                    if (scalor[i] != 0)
                    {
                        scaledArray[i] = array[i] / scalor[i];
                    }
                    else
                    {
                        // basic fix to prevent division by zero
                        scalor[i] = 0.00001;
                        scaledArray[i] = array[i] / scalor[i];

                    }
                }
            }

            return scaledArray;
        }

        public static double[,] Diag(this double[,] matrix, double[] diagVals)
        {
            int rows = matrix.GetLength(0);
            int cols = matrix.GetLength(1);
            // the matrix has to be scare
            if (rows == cols)
            {
                double[,] diagMatrix = new double[rows, cols];
                int k = 0;
                for (int i = 0; i < rows; i++)
                    for (int j = 0; j < cols; j++)
                    {
                        if (i == j)
                        {
                            diagMatrix[i, j] = diagVals[k];
                            k++;
                        }
                        else
                        {
                            diagMatrix[i, j] = 0;
                        }
                    }
                return diagMatrix;
            }
            else
            {
                Console.WriteLine("Diag should be used on square matrix only.");
                return null;
            }


        }
    }

    #endregion
    #region FOUNDATION
    //internal interface IInterpolate
    //{
    //    double? Interpolate(double p);
    //}

    internal abstract class Interpolation //: IInterpolate
    {

        public Interpolation(double[] _x, double[] _y)
        {
            int xLength = _x.Length;
            if (xLength == _y.Length && xLength > 1 && _x.Distinct().Count() == xLength)
            {
                x = _x;
                y = _y;
            }
        }

        // cubic spline relies on the abscissa values to be sorted
        public Interpolation(double[] _x, double[] _y, bool checkSorted = true)
        {
            int xLength = _x.Length;
            if (checkSorted)
            {
                if (xLength == _y.Length && xLength > 1 && _x.Distinct().Count() == xLength && Enumerable.SequenceEqual(_x.SortedList(), _x.ToList()))
                {
                    x = _x;
                    y = _y;
                }
            }
            else
            {
                if (xLength == _y.Length && xLength > 1 && _x.Distinct().Count() == xLength)
                {
                    x = _x;
                    y = _y;
                }
            }
        }

        public double[] X
        {
            get
            {
                return x;
            }
        }

        public double[] Y
        {
            get
            {
                return y;
            }
        }

        public abstract double? Interpolate(double p);

        private double[] x;
        private double[] y;
    }

    #endregion

    #region LINEAR
    internal sealed class LinearInterpolation : Interpolation
    {

        public LinearInterpolation(double[] _x, double[] _y) : base(_x, _y)
        {
            len = X.Length;
            if (len > 1)
            {
                Console.WriteLine("Successfully set abscissa and ordinate.");
                baseset = true;
                // make a copy of X as a list for later use
                lX = X.ToList();
            }
            else
                Console.WriteLine("Ensure x and y are the same length and have at least 2 elements. All x values must be unique.");
        }

        public override double? Interpolate(double p)
        {
            if (baseset)
            {
                double? result = null;
                double Rx;

                try
                {
                    // point p may be outside abscissa's range
                    // if it is, we return null
                    Rx = X.First(s => s >= p);
                }
                catch (ArgumentNullException)
                {
                    return null;
                }

                // at this stage we know that Rx contains a valid value
                // find the index of the value close to the point required to be interpolated for
                int i = lX.IndexOf(Rx);

                // provide for index not found and lower and upper tabulated bounds
                if (i == -1)
                    return null;

                if (i == len - 1 && X[i] == p)
                    return Y[len - 1];

                if (i == 0)
                    return Y[0];

                // linearly interpolate between two adjacent points
                double h = (X[i] - X[i - 1]);
                double A = (X[i] - p) / h;
                double B = (p - X[i - 1]) / h;

                result = Y[i - 1] * A + Y[i] * B;

                return result;

            }
            else
            {
                return null;
            }
        }

        private bool baseset = false;
        private int len;
        private List<double> lX;
    }

    #endregion

    #region  CUBIC_SPLINE
    internal sealed class CubicSplineInterpolation : Interpolation
    {
        public CubicSplineInterpolation(double[] _x, double[] _y) : base(_x, _y, true)
        {
            len = X.Length;
            if (len > 1)
            {

                Console.WriteLine("Successfully set abscissa and ordinate.");
                baseset = true;
                // make a copy of X as a list for later use
                lX = X.ToList();
            }
            else
                Console.WriteLine("Ensure x and y are the same length and have at least 2 elements. All x values must be unique. X-values must be in ascending order.");
        }

        public override double? Interpolate(double p)
        {
            if (baseset)
            {
                double? result = 0;
                int N = len - 1;

                double[] h = X.Diff();
                double[] D = Y.Diff().Scale(h, false);
                double[] s = Enumerable.Repeat(3.0, D.Length).ToArray();
                double[] dD3 = D.Diff().Scale(s);
                double[] a = Y;

                // generate tridiagonal system
                double[,] H = new double[N - 1, N - 1];
                double[] diagVals = new double[N - 1];
                for (int i = 1; i < N; i++)
                {
                    diagVals[i - 1] = 2 * (h[i - 1] + h[i]);
                }

                H = H.Diag(diagVals);

                // H can be null if non-square matrix is passed
                if (H != null)
                {
                    for (int i = 0; i < N - 2; i++)
                    {
                        H[i, i + 1] = h[i + 1];
                        H[i + 1, i] = h[i + 1];
                    }

                    double[] c = new double[N + 2];
                    c = Enumerable.Repeat(0.0, N + 1).ToArray();

                    // solve tridiagonal matrix
                    TridiagonalMatrix.Tridiagonal le = new TridiagonalMatrix.Tridiagonal();
                    double[] solution = le.Solve(H, dD3);

                    for (int i = 1; i < N; i++)
                    {
                        c[i] = solution[i - 1];
                    }

                    double[] b = new double[N];
                    double[] d = new double[N];
                    for (int i = 0; i < N; i++)
                    {
                        b[i] = D[i] - (h[i] * (c[i + 1] + 2.0 * c[i])) / 3.0;
                        d[i] = (c[i + 1] - c[i]) / (3.0 * h[i]);
                    }

                    double Rx;

                    try
                    {
                        // point p may be outside abscissa's range
                        // if it is, we return null
                        Rx = X.First(m => m >= p);
                    }
                    catch
                    {
                        return null;
                    }

                    // at this stage we know that Rx contains a valid value
                    // find the index of the value close to the point required to be interpolated for
                    int iRx = lX.IndexOf(Rx);

                    if (iRx == -1)
                        return null;

                    if (iRx == len - 1 && X[iRx] == p)
                        return Y[len - 1];

                    if (iRx == 0)
                        return Y[0];

                    iRx = lX.IndexOf(Rx) - 1;
                    Rx = p - X[iRx];
                    result = a[iRx] + Rx * (b[iRx] + Rx * (c[iRx] + Rx * d[iRx]));

                    return result;
                }
                else
                {
                    return null;
                }


            }
            else
                return null;
        }

        private bool baseset = false;
        private int len;
        private List<double> lX;
    }
    #endregion
    class EntryPoint
    {
        static void Main(string[] args)
        {
            // code relies on abscissa values to be sorted
            // there is a check for this condition, but no fix
            // f(x) = 1/(1+x^2)*sin(x)


            string[] MeasurementSpectrumData;   // 측정 스펙트럼 데이터 저장할 배열. (한 줄씩 저장)
            string[] SingleLineData;            // 한 줄의 스펙트럼 데이터를 임시로 저장할 배열.

            // "SiO2 2nm_on_Si.dat" 파일 읽기. (한 줄씩)
            MeasurementSpectrumData = File.ReadAllLines("Si_nm.txt");
            int LoopNum = MeasurementSpectrumData.Length;

            // wavelength : 350 ~ 1000(nm)인 측정 스펙트럼 데이터를 담을 리스트 선언.
            //List<double> wavelength = new List<double>();   // 파장 데이터 리스트.
            //List<double> n = new List<double>();    // n 데이터 리스트.
            //List<double> k = new List<double>();   // k 데이터 리스트.
            double[] wavelength = new double[LoopNum-1];
            double[] n = new double[LoopNum-1];
            double[] k = new double[LoopNum-1];

            // 데이터의 첫번째 줄은 column 명이다.
            // 이를 제외하기 위해 반복문을 1부터 시작한다.
            int StartIndex = 1;
            
            for (int i = StartIndex; i < LoopNum; i++)
            {
                // tsv 형식의 데이터를 SingleLineData에 저장한다.
                SingleLineData = MeasurementSpectrumData[i].Split((char)0x09);  // 0x09 : 수평 탭.

                // 각 컬럼에 해당하는 데이터를 저장한다.
                //wavelength.Add(Double.Parse(SingleLineData[0]));
                //n.Add(Double.Parse(SingleLineData[1]));
                //k.Add(Double.Parse(SingleLineData[2]));
                wavelength[i-1] = double.Parse(SingleLineData[0]);
                n[i-1] = double.Parse(SingleLineData[1]);
                k[i-1] = double.Parse(SingleLineData[2]);
            }

            for (int i = 0; i < LoopNum-1; i++)
            {
                Console.WriteLine($"{wavelength[i]}          {n[i]}          {k[i]}");
            }
            /*double[] ex = new double[wavelength.Count];
            double[] en = new double[wavelength.Count];
            //double[] p = new double[wavelength.Count];
            double[] p = { 350, 355, 360, 400, 450, 500, 600, 650, 700, 750, 800, 850 };

            for (int i = 0; i < wavelength.Count; i++)
            {
                ex[i] = wavelength[i];
                en[i] = n[i];
            }

            LinearInterpolation LI = new LinearInterpolation(ex, en);
            CubicSplineInterpolation CS = new CubicSplineInterpolation(ex, en);

            Console.WriteLine("Linear Interpolation:");
            foreach (double pp in p)
            {
                Console.WriteLine($"{pp} \t:\t"
                + LI.Interpolate(pp).ToString());

            }
            Console.WriteLine();
            Console.WriteLine();

            Console.WriteLine("Cubic Spline Interpolation:");
            foreach (double pp in p)
            {

                Console.WriteLine($"{pp} \t:\t"
                    + CS.Interpolate(pp).ToString());
            }
            Console.ReadLine();*/
        }
    }
}