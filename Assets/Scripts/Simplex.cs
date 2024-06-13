namespace DTSimulation
{
    public class Simplex
    {
        // opposite neighbors of triangle relative to each vertex
        // ie: for vertex A, the triangle against the opposite edge of that vertex is its neighbor
        //     X
        //    / \
        //   /   \
        //  /  ^  \
        // X---|---X
        //  \  |  / 
        //   \ | /
        //     A
        //
        public Simplex[] neighbors;
        public int[] vertices;
        public int sum;
        public int label;
        public bool flag;

        public Simplex(int[] a, int N, int pointer_number)
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
            label = pointer_number;
            flag = false;
        }
    }
}
