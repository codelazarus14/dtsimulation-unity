using UnityEngine;

namespace DTSimulation
{
    public class Embed
    {
        public float[] X { get; private set; }
        public float[] Y { get; private set; }
        public int[] Mark { get; private set; }

        private DT dt;
        private int n;

        public Embed(DT dt)
        {
            this.dt = dt;
            // simplex_number=pointer_number after tidy()
            n = dt.node_number;
            X = new float[n];
            Y = new float[n];
        }

        public void ComputeEmbedding()
        {
            int coord;
            FormMatrix();
            coord = 0;
            CG_Solver(coord);
            coord = 1;
            CG_Solver(coord);
        }

        private void SparseMult(float[] v, float[] p)
        {
            for (int i = 0; i < n; i++)
            {
                if (Mark[i] == 0)
                {
                    p[i] = 0f;
                    for (int j = dt.nstart[i]; j < dt.nstart[i + 1]; j++)
                    {
                        if (Mark[dt.ncol[j]] == 0)
                            p[i] = p[i] + dt.nlap[j] * v[dt.ncol[j]];
                    }
                }
            }
        }

        private void FormMatrix()
        {
            Mark = new int[n];

            for (int i = 0; i < n; i++)
                Mark[i] = 0;

            Mark[0] = dt.boundary_length + 1;
            for (int i = 0; i < dt.boundary_length; i++)
                Mark[dt.boundary[i]] = i + 1;
        }

        private float Dot(float[] a, float[] b)
        {
            float d = 0f;
            for (int i = 0; i < n; i++)
            {
                if (Mark[i] == 0)
                    d += a[i] * b[i];
            }

            return d;
        }

        // remove one vertex from sphere. This makes disk.
        // Pin neighbour triangles around circle and solve for
        // interior points placing each at center of mass of its
        // neighbours
        private void CG_Solver(int coord)
        {
            float[] R_0 = new float[n];
            float[] R_1 = new float[n];
            float[] P_0 = new float[n];
            float[] P_1 = new float[n];
            float[] SS = new float[n];
            float[] SOL = new float[n];
            float[] SOL_1 = new float[n];
            float[] B = new float[n];

            float alpha, beta, rrtmp, rrtmp2, resid, psdot;
            int count = 0, itmp, done;
            float CG_RESIDUAL = 0.0000001f;

            // pin neighbours of simplex 0 to unit circle. Yields RHS/source for
            // laplacian problem

            for (int i = 0; i < n; i++)
                B[i] = 0f;

            for (int i = 0; i < n; i++)
            {
                if (Mark[i] == 0)
                {
                    // bulk pt
                    for (int j = dt.nstart[i]; j < dt.nstart[i + 1]; j++)
                    {
                        itmp = dt.ncol[j];
                        if (Mark[itmp] != 0)
                        {
                            // connects to boundary pt
                            if (coord == 0)
                                B[i] = B[i] + Mathf.Cos((2 * Mathf.PI / dt.boundary_length) * Mark[itmp]);
                            else if (coord == 1)
                                B[i] = B[i] + Mathf.Sin((2 * Mathf.PI / dt.boundary_length) * Mark[itmp]);
                        }
                    }
                }
            }

            for (int i = 0; i < n; i++)
            {
                SOL[i] = 0f;
                P_0[i] = B[i];
                R_0[i] = B[i];
            }

            do
            {
                SparseMult(P_0, SS);

                rrtmp = Dot(R_0, R_0);
                psdot = Dot(SS, P_0);

                alpha = rrtmp / psdot;

                for (int i = 1; i < n; i++)
                    R_1[i] = R_0[i] - alpha * SS[i];

                for (int i = 1; i < n; i++)
                    SOL_1[i] = SOL[i] + alpha * P_0[i];

                rrtmp2 = Dot(R_1, R_1);

                beta = rrtmp2 / rrtmp;

                for (int i = 1; i < n; i++)
                    P_1[i] = R_1[i] + beta * P_0[i];
                resid = Mathf.Sqrt(rrtmp2 / n);


                for (int i = 1; i < n; i++)
                {
                    R_0[i] = R_1[i];
                    P_0[i] = P_1[i];
                    SOL[i] = SOL_1[i];
                }

                count++;
            } while ((resid > CG_RESIDUAL) && (count < 1000));


            for (int i = 1; i < n; i++)
            {
                if (coord == 0)
                    X[i] = SOL[i];
                else if (coord == 1)
                    Y[i] = SOL[i];
            }

            for (int j = 0; j < dt.boundary_length; j++)
            {
                if (coord == 0)
                    X[dt.boundary[j]] = Mathf.Cos((2 * Mathf.PI / dt.boundary_length) * (j + 1));
                else if (coord == 1)
                    Y[dt.boundary[j]] = Mathf.Sin((2 * Mathf.PI / dt.boundary_length) * (j + 1));
            }
        }
    }
}
