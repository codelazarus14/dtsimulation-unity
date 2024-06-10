namespace DTSimulation
{
    public class Simplex
    {
        public Simplex[] neighbors;
        public int[] vertices;
        public int sum;
        public int label;
        public bool flag;

        public Simplex(int[] a, int N)
        {
            int i, count;
            int DPLUS = N + 1;

            vertices = new int[DPLUS];
            neighbors = new Simplex[DPLUS];
            count = 0;

            for (i = 0; i < DPLUS; i++)
            {
                vertices[i] = a[i];
                neighbors[i] = null;
                count += a[i];
            }

            sum = count;
            flag = false;
        }
    }
}
