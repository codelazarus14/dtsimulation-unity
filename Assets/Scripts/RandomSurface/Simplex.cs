namespace DTSimulation.RandomSurface
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
        // array of vertex labels
        public int[] vertices;
        public int sum;
        public int label;
        public bool flag;

        public Simplex(int[] a, int N, int pointer_number)
        {
            vertices = new int[N];
            neighbors = new Simplex[N];
            sum = 0;

            for (int i = 0; i < N; i++)
            {
                vertices[i] = a[i];
                neighbors[i] = null;
                sum += a[i];
            }

            label = pointer_number;
            flag = false;
        }
    }
}
