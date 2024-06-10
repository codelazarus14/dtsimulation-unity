using UnityEngine;

public class DTSimulation : MonoBehaviour
{
    // constants
    public const int DISKVOL = 472;
    public const int DISKMINVOL = DISKVOL - 30;
    public const int DISKMAXVOL = DISKVOL + 30;
    public const int MARKEDQ = 282; //FRAC*(DISKVOL+2)/(2.0-FRAC)

    public const int VOL = DISKVOL * MARKEDQ;
    public const int BIGVOL = 4 * VOL;

    [SerializeField]
    private int dimensions;

    [SerializeField]
    private Simplex[] simplexes;

    private void Start() => Main();

    void Main()
    {
        // abstract graph of points
        // to some coordinate space = enforce smoothness = flatness thru an Action:
        // - maximize dot product of 3d triangle normals? parallel/coplanar
        // - embedded curvature
        // take curr action value, imagine doing it, compute new action resulting, if < old? accept, otherwise accept stochastically... exponential
        // moves towards lower energies with some fluctuations towards higher energy - metropolis algo

        // pick moves at random, link flips or moving vertices
        // generally move towards lower energy w some fluctuations towards higher (depending on temperature)

        int iter, sweep;

        Thermalize();
    }

    private void Thermalize()
    {
        int therm;

        simplexes = new Simplex[BIGVOL];


    }

    private void ReadFile()
    {

    }
}