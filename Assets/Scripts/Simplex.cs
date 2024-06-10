public class Simplex
{
    public Simplex[] neighbors;
    public int[] vertex;
    public int sum;
    public int label;
    public bool flag;

    public Simplex(int[] a, int N)
    {
        int i, count;
        int DPLUS = N + 1;

        vertex = new int[DPLUS];
        neighbors = new Simplex[DPLUS];
        count = 0;

        for (i = 0; i < DPLUS; i++)
        {
            vertex[i] = a[i];
            neighbors[i] = null;
            count += a[i];
        }

        sum = count;
        flag = false;
    }
}