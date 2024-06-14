using System;
using System.Collections.Generic;
using UnityEngine;
using Random = UnityEngine.Random;

namespace DTSimulation
{
    public class DT : MonoBehaviour
    {
        // all comments using /**/ are from DT.java

        // TODO: fix style/naming conventions
        // TODO: fix messy/redundant loops and reduce work
        // TODO: remove unused, add cached values, make LOGICS and ints = 0 or 1 become bools (usually "int done")

        /* some global constants for run */

        private float kappa_0 = 0f, kappa_d = 0f, kappa_0b = 0f;
        private const float BETA = 3f;
        private const float ALPHA = 0.25f;
        private const int DISKVOL = 472;
        private const int DISKMINVOL = DISKVOL - 30;
        private const int DISKMAXVOL = DISKVOL + 30;
        private const float FRAC = 0.75f;
        private const float TC = 7f;
        private const int MARKEDQ = 282; //FRAC*(DISKVOL+2)/(2.0-FRAC)

        private const int D = 2;
        private const int DPLUS = D + 1;
        private const int DPLUSPLUS = D + 2;
        private const int VOL = DISKVOL + MARKEDQ;
        private const int BIGVOL = 4 * VOL;
        private const int MAXVOL = VOL + 30;
        private const int MINVOL = DPLUSPLUS;
        private const int NUMROW = VOL / 2;
        private const int NONZEROS = NUMROW * VOL;

        private const int THERMALISE = 1;
        private const int SWEEPS = 100;
        private const int TUNE_COUPLING = 100;
        private const int SEED = 1;
        private const int GAP = 10;
        private const int DV = 2;

        /* simple global pointers and counters */

        private Node stack_head = null;
        private int simplex_number = 0;
        private int pointer_number = 0;
        private int boundary_length;
        private int node_number = 0;

        private int[] boundary;

        private int[] nstart, ncol;
        private float[] nlap;

        /* data structures */

        private Simplex[] simplex_point;
        private int stack_count = 0;

        /* measurements etc */

        private int number_measure = 0, max_point;
        private int[] legal_subsimplex, try_subsimplex, manifold_subsimplex, go_subsimplex;
        private int vol_monitor = 0, num_monitor = 0, b_monitor = 0, growing_vol;
        private float manifold_check;

        private bool grow;

        /* simple observables -- number of nodes, simplices and size */
        /* and variable versions of couplings */

        private float real_simplex = 0f, real_node = 0f;

        /* routines checks manifold constraint */
        /* this entails examining neighbour simplices for an occurrence */
        /* of the 'opposing' vertex in a simplex which also contains the other*/
        /* new common subsimplex vertices */
        /* this is equivalent to requiring that the new common subsimplex */
        /* is not present anywhere else in the triangulation */
        /* examine all neighbours gotten by moving out on faces */
        /* which contain the other new common subsimplex vertices */

        //
        // my fields/deviations from above
        //

        [SerializeField]
        private TextAsset config;


        private void Start() => Main();

        /* driver for simplex Monte Carlo */
        private void Main()
        {
            // abstract graph of points
            // to some coordinate space = enforce smoothness = flatness thru an Action:
            // - maximize dot product of 3d triangle normals? parallel/coplanar
            // - embedded curvature
            // take curr action value, imagine doing it, compute new action resulting, if < old? accept, otherwise accept stochastically... exponential
            // moves towards lower energies with some fluctuations towards higher energy - metropolis algo

            // pick moves at random, link flips or moving vertices
            // generally move towards lower energy w some fluctuations towards higher (depending on temperature)

            // TODO: make thermalize a coroutine and set value back up to 500, sweeps to 20000
            Thermalize();
            /* sweep lattice outputting measurements every so often */

            Debug.Log("Starting measurements");
            for (int sweep = 1; sweep < SWEEPS; sweep++)
            {
                Debug.Log($"sweep: {sweep}");
                for (int iter = 0; iter < VOL; iter++)
                    TrialChange();
                Tidy();

                bool done = false;
                int v;
                Debug.Log($"Simplex number: {simplex_number}");
                Simplex p = null;
                int[] num = new int[1];
                int[] nn = new int[VOL];

                for (int i = 0; i < simplex_number; i++)
                {
                    if (done) break;
                    for (int j = 0; j < DPLUS; j++)
                    {
                        if (simplex_point[i].vertices[j] == 0)
                        {
                            p = simplex_point[i];
                            done = true;
                            break;
                        }
                    }
                }
                v = 0;

                GetNN(p, v, nn, num);

                boundary_length = num[0];

                num_monitor++;
                vol_monitor += simplex_number;
                b_monitor += boundary_length;

                if (sweep % 100 == 0)
                {
                    ShiftCoupling();
                    Debug.Log($"Boundaryyy length: {boundary_length}");
                    Debug.Log($"Total bulkyyy node: {node_number - boundary_length - 1}");
                }

                if (sweep % GAP == 0)
                    Measure();
            }

            /* finally dump final configuration and print some results */
            PrintOut();
        }

        private void Thermalize()
        {
            simplex_point = new Simplex[BIGVOL];
            ReadFile();
            //InitialConfig();

            nlap = new float[NONZEROS];
            ncol = new int[NONZEROS];
            nstart = new int[VOL];
            boundary = new int[VOL];

            //PrintConfig();

            legal_subsimplex = new int[DPLUS];
            try_subsimplex = new int[DPLUS];
            manifold_subsimplex = new int[DPLUS];
            go_subsimplex = new int[DPLUS];

            /* build lattice */


            //grow = LOGIC.NO;
            //if (grow == LOGIC.YES)
            //{
            //    while (growing_vol < VOL)
            //    {
            //        //System.out.println(growing_vol);
            //        trial_change();
            //        growing_vol += D;
            //        //PrintConfig();
            //    }
            //}

            //Tidy();

            grow = false;

            /* thermalise and output info on run */

            Header();

            Debug.Log("Thermalizing lattice");

            for (int therm = 1; therm < THERMALISE; therm++)
            {
                for (int iter = 0; iter < VOL; iter++)
                    TrialChange();
                Tidy();

                // find pointer to node 0
                int done = 0, v;
                Simplex p = null;
                int[] nn = null;
                int[] num = null;
                num = new int[1];
                nn = new int[VOL];


                for (int i = 0; i < simplex_number; i++)
                {
                    if (done == 1) break;
                    for (int j = 0; j < (D + 1); j++)
                    {
                        if (simplex_point[i].vertices[j] == 0) { p = simplex_point[i]; done = 1; break; }
                    }
                }
                v = 0;

                GetNN(p, v, nn, num);

                num_monitor++;
                vol_monitor += simplex_number;
                b_monitor += num[0];

                string message = "";
                if (therm % TUNE_COUPLING == 0)
                {
                    message += "V: " + (simplex_number - num[0]);
                    message += "L: " + num[0];
                    message += "N0B: " + (node_number - boundary_length - 1);
                    Debug.Log(message);

                    ShiftCoupling();
                }
            }

            Init();
        }

        private void PrintConfig()
        {
            string message = "Configuration\n";
            message += "Number of simplices and nodes " + simplex_number + ", " + node_number + "\n";
            for (int i = 0; i < pointer_number; i++)
            {
                message += "Simplex " + i + "\n";
                if (simplex_point[i] != null)
                {
                    for (int j = 0; j < DPLUS; j++)
                    {
                        message += " - vertex " + simplex_point[i].vertices[j] + "\n";
                        message += " - neighbor " + simplex_point[i].neighbors[j].label + "\n";
                    }
                }
            }
            Debug.Log(message);
        }

        private bool AllowedMove(Simplex p, int sub, int[] a)
        {
            Simplex[] examine = new Simplex[VOL];
            Simplex[] array1 = new Simplex[VOL];
            Simplex[] array2 = new Simplex[VOL];
            Simplex[] dum = new Simplex[DPLUS];
            int[] b = new int[DPLUS];

            int j, number1, number2, search;

            bool good = true;

            /* need to teminate when either 1. see opposing vertex a[DPLUS] in some */
            /* simplex or 2. have examined all simplices containing the new common */
            /* vertices derived from original simplex */
            /* so need to flag simplices that have been examined */

            if (sub == D)
                return good;
            if (sub == 0)
                return good;

            j = 0;
            for (int i = 0; i < DPLUS; i++)
            {
                if (i > sub)
                {
                    b[j] = a[i];
                    j++;
                }
            }

            array1[0] = p;
            number1 = 1;
            examine[0] = p;
            search = 1;
            p.flag = true;

            /* loop while new neighbour simplices to examine */

            while (number1 > 0)
            {
                number2 = 0;

                for (int i = 0; i < number1; i++)
                {
                    /* examine to see if contains 'opposing vertex */

                    for (j = 0; j < DPLUS; j++)
                    {
                        if (array1[i].vertices[j] == a[DPLUS])
                        {
                            good = false;
                            break;
                        }
                    }

                    if (!good) break;

                    /* find simplices which also in common with this subsimplex */
                    /* and neighbour to simplex array1[i] */

                    CommonSimplex(array1[i], b, D - sub, dum);

                    /* check to see whether any of these have been seen before */
                    /* if not then add to array2 and flag seen. Also make note in examine */

                    for (j = 0; j < sub + 1; j++)
                    {
                        if (!dum[j].flag)
                        {
                            array2[number2] = dum[j];
                            number2++;
                            dum[j].flag = true;

                            examine[search] = dum[j];
                            search++;
                        }
                    }
                }

                if (!good) break;

                for (int i = 0; i < number2; i++)
                    array1[i] = array2[i];

                number1 = number2;
            }

            /* depending on value of good we have either found the opposing vertex or */
            /* examined all the associated simplices */

            /* first set all local flags back to zero */

            for (int i = 0; i < search; i++)
                examine[i].flag = false;

            manifold_check += (float)search / VOL;
            return good;
        }

        /* routine takes n element vector a[] and returns n-1 element in b[] */
        /* thus n combinations possible selected according to count */
        /* count is the index 'left out' in getting n-1 from n */
        /* vertex left out of final vector b[] is returned in b[n-1] */
        private void Combo(int[] a, int[] b, int n, int count)
        {
            int add = 0;
            for (int i = 0; i < n; i++)
            {
                if (count != i)
                {
                    b[add] = a[i];
                    add++;
                }
                else b[n - 1] = a[i];
            }
            return;
        }

        /* routine takes a named subsimplex a[0]...a[n-1] and searches for it */
        /* in a simplex pointed at by p. */
        /* returns a vector containing the local indices of pointers to */
        /* neighbouring simplices */
        /* which share a face also encompassing this subsimplex */
        private void CommonSimplex(Simplex p, int[] a, int n, Simplex[] face)
        {
            int j = 0;
            int[] b = new int[DPLUS];
            int[] mask = new int[DPLUS];
            bool[] found = new bool[DPLUS];

            for (int i = 0; i < n; i++)
                found[i] = false;

            /* find positions/local indices of subsimplex in this simplex */

            for (int i = 0; i < n; i++)
            {
                for (j = 0; j < DPLUS; j++)
                {
                    if (p.vertices[j] == a[i])
                    {
                        b[i] = j;
                        found[i] = true;
                        break;
                    }
                }
            }

            for (int i = 0; i < n; i++)
            {
                if (!found[i])
                {
                    Debug.LogError("Error in CommonSimplex");
                    Application.Quit(); // replacement for System.exit(1);
                }
            }

            for (j = 0; j < DPLUS; j++)
                mask[j] = 0;

            for (j = 0; j < n; j++)
                mask[b[j]] = 1;

            j = 0;
            for (int i = 0; i < DPLUS; i++)
            {
                if (mask[i] == 0)
                {
                    face[j] = p.neighbors[i];
                    j++;
                }
            }

            return;
        }

        /* routine takes pointer to simplex and a vertex and returns local */
        /* index to conjugate face */
        private int FindFace(Simplex p, int a)
        {
            for (int i = 0; i < DPLUS; i++)
                if (p.vertices[i] == a) return i;

            Debug.LogError("Error in FindFace");
            Application.Quit(); // replacement for System.exit(1);

            /* return dummy if get to here error */
            return VOL;
        }

        // TODO: about comment below, refactor FindSimplices to just return the value instead?
        // do we even need FindOrder? it just calls FindSimplices and returns the value lmao

        // wrap num in single element array to pass by ref !    <-- from java
        private int FindOrder(int a, Simplex pp)
        {
            Simplex[] list = new Simplex[MAXVOL];
            int[] dum = new int[DPLUS];
            int[] num = new int[1];
            /* finds num of d-simps sharing vertex a */
            dum[0] = a;
            int dummy = 1;
            FindSimplices(pp, dum, dummy, list, num);
            return num[0];
        }

        /* finds addresses of all simplices which share a given subsimplex */
        private void FindSimplices(Simplex p, int[] a, int sub, Simplex[] s_near, int[] num)
        {
            // TODO: a lot of this code seems to be copied from AllowedMove
            int num1, num2;
            Simplex[] array1 = new Simplex[MAXVOL];
            Simplex[] array2 = new Simplex[MAXVOL];
            Simplex[] near = new Simplex[DPLUS];

            int number;

            array1[0] = p;
            num1 = 1;
            s_near[0] = p;
            number = 1;
            s_near[0].flag = true;

            while (num1 > 0)
            {
                num2 = 0;

                for (int i = 0; i < num1; i++)
                {
                    CommonSimplex(array1[i], a, sub, near);

                    for (int j = 0; j < DPLUS - sub; j++)
                    {
                        if (!near[j].flag)
                        {
                            s_near[number] = near[j];
                            array2[num2] = near[j];
                            near[j].flag = true;
                            number++;
                            num2++;
                        }
                    }
                }

                for (int i = 0; i < num2; i++)
                    array1[i] = array2[i];

                num1 = num2;
            }

            for (int i = 0; i < number; i++)
                s_near[i].flag = false;

            num[0] = number;
            return;
        }

        /* takes simplex and order of subsimplex returns logical flag if legal move */
        /* legal move is equivalent to being exactly d+1-i associated simplices. */
        /* these subsimplex vertices occupy first i+1 entries in a[] */
        /* other d-i 'external' vertices of original simplex occupy end of a[] */
        /* also returns pointers to these d+1-i simplices */
        /* and opposing vertex is placed at end of a[] */
        private bool GoodSubsimplex(Simplex p, int sub, int[] a, Simplex[] isimplex)
        {
            int add, temp, opposing;
            bool seen_already;
            int[] aind = new int[DPLUS];

            /* test whether subsimplex is simplex itself i.e node insertion move */

            if (sub == D)
            {
                for (int i = 0; i < DPLUS; i++)
                    a[i] = p.vertices[i];

                /* opposing vertex in this case is 'new' one obtained off stack */

                if (stack_head == null)
                    a[DPLUS] = node_number;
                else
                    a[DPLUS] = stack_head.name;

                isimplex[DPLUS] = p;
                return true;
            }

            /* otherwise generate subsimplex at random placing its indices in aind */

            add = 0;
            while (add < sub + 1)
            {
                temp = Random.Range(0, DPLUS);
                /* generate random index */

                /* scan existing ones to see if already produced */

                seen_already = false;
                for (int i = 0; i < add; i++)
                {
                    if (temp == aind[i])
                    {
                        seen_already = true;
                        break;
                    }
                }

                if (!seen_already)
                {
                    aind[add] = temp;
                    add++;
                }
            }

            /* now create array of indices to d-i remaining vertices */

            temp = add;
            for (int i = 0; i < DPLUS; i++)
            {
                seen_already = false;
                for (add = 0; add < sub + 1; add++)
                {
                    if (i == aind[add])
                    {
                        seen_already = true;
                        break;
                    }
                }

                if (!seen_already)
                {
                    aind[temp] = i;
                    temp++;
                }
            }

            if (temp != DPLUS)
            {
                Debug.LogError("Error in GoodSubsimplex");
                Application.Quit(); // replacement for System.exit(1);
            }


            /* now loop over all possible faces constucted to include this subsimplex */
            /* by selecting d-i-1 out of the d-i remaining indices */

            for (int i = 0; i < DPLUS; i++)
                a[i] = p.vertices[aind[i]];

            isimplex[sub + 1] = p.neighbors[aind[sub + 1]];
            opposing = isimplex[sub + 1].sum - SumFace(p, a[sub + 1]);

            /* protect following loop if sub=D-1 */

            if (sub < (D - 1))
            {
                for (int i = sub + 2; i < DPLUS; i++)
                {
                    isimplex[i] = p.neighbors[aind[i]];
                    if (isimplex[i].sum - SumFace(p, a[i]) != opposing)
                        return false;
                }
            }

            a[DPLUS] = opposing;
            isimplex[DPLUS] = p;
            return true;
        }

        private void ReadFile()
        {
            // TODO: store state in some kind of struct that can be serialized/inspected
            // the temp field, for example, is literally just throwing away some of the Scanner-esque outputs
            // bc those fields are only edited at runtime, not loaded in
            // basically this is all based on a fragile text-based structure that has to be 1-1 with the whitespaced config file
            // when using a structured JSON type thing would save a lot of headaches and look better

            // also TODO: maybe create a textmesh to store these outputs like an in-game console cuz the debug one is not great

            // load config
            string[] configLines = config.text.Split(new char[0], StringSplitOptions.RemoveEmptyEntries);
            // using a string stack to imitate use of the Scanner in the java version
            Stack<string> configStack = new Stack<string>();
            for (int i = configLines.Length - 1; i >= 0; i--)
            {
                configStack.Push(configLines[i]);
            }

            int stackCount, s_number, dummy, l_number, c, e1, e2, simp;

            int[] dum = new int[DPLUS];
            int[,] dum2;

            float k0, b;
            float temp;

            // read config values
            try
            {
                s_number = int.Parse(configStack.Pop());
                dum2 = new int[s_number, DPLUS];
                node_number = int.Parse(configStack.Pop());
                stackCount = int.Parse(configStack.Pop());
                // TODO: why do these exist + debug cleanup
                temp = float.Parse(configStack.Pop());
                temp = float.Parse(configStack.Pop());
                temp = float.Parse(configStack.Pop());

                Debug.Log($"{s_number}, {node_number}, {stackCount}");

                simplex_number = 0;
                stack_count = 0;
                pointer_number = 0;

                Debug.Log("Reading in existing configuration");

                Debug.Log($"Configuration has volume: {s_number}");
                Debug.Log($"Node number: {node_number}");
                Debug.Log($"k0b coupling, k2 coupling and curvaturesq coupling: {kappa_0b}, {kappa_d}, {BETA}");

                for (int i = 0; i < stackCount; i++)
                {
                    // add used vertex labels to stack
                    dummy = int.Parse(configStack.Pop());
                    Push(dummy);
                }

                // set up simplices
                for (int i = 0; i < s_number; i++)
                {
                    for (int j = 0; j < DPLUS; j++)
                    {
                        dum[j] = int.Parse(configStack.Pop());
                        dum2[i, j] = int.Parse(configStack.Pop());
                    }

                    simplex_point[i] = new Simplex(dum, D, pointer_number);
                    simplex_number++;
                    pointer_number++;
                }

                // set neighbors
                for (int i = 0; i < s_number; i++)
                {
                    for (int j = 0; j < DPLUS; j++)
                    {
                        simplex_point[i].neighbors[j] = simplex_point[dum2[i, j]];
                    }
                }

                // print out results
                string message = "";
                for (int i = 0; i < s_number; i++)
                {
                    for (int j = 0; j < DPLUS; j++)
                    {
                        message += simplex_point[i].vertices[j] + "\t" + simplex_point[i].neighbors[j].label + "\t";
                    }
                    message += "\n";
                }
                Debug.Log(message);

                growing_vol = simplex_number;

                Debug.Log("Have read data successfully");
            }
            catch (FormatException e)
            {
                Debug.LogError($"Error parsing config file! {e.Message}");
            }
        }

        // TODO: rename "GetNearestNeighbor"

        // hardwired for D=2 right now ..
        private void GetNN(Simplex p, int v, int[] nn, int[] vnum)
        {
            Simplex[] list = new Simplex[BIGVOL];
            int[] dum = new int[DPLUS];
            int[] seen = new int[VOL];
            int[] num = new int[1];
            int[] v1 = new int[VOL];
            int[] v2 = new int[VOL];
            int k, currentpt, nextpt, working, index;

            dum[0] = v;
            int dummy = 1;
            FindSimplices(p, dum, dummy, list, num);

            for (int i = 0; i < VOL; i++)
                seen[i] = 0;

            k = 0;
            for (int i = 0; i < num[0]; i++)
            {
                index = 0;

                for (int j = 0; j < DPLUS; j++)
                {
                    if (list[i].vertices[j] == v)
                        index = j;
                }

                v1[k] = list[i].vertices[(index + 1) % DPLUS];
                v2[k] = list[i].vertices[(index + 2) % DPLUS];
                k++;
            }

            nn[0] = v1[0];
            seen[0] = 1;
            k = 1;

            // look for v1[0] in rest of v1/v2 arrays
            currentpt = v1[0];
            do
            {
                nextpt = 0;
                for (int i = 0; i < num[0]; i++)
                {
                    if (seen[i] == 0)
                    {
                        if (v1[i] == currentpt)
                        {
                            nn[k] = v2[i];
                            nextpt = v2[i];
                            k++;
                            seen[i] = 1;
                        }
                        if (v2[i] == currentpt)
                        {
                            nn[k] = v1[i];
                            nextpt = v1[i];
                            k++;
                            seen[i] = 1;
                        }
                    }
                }
                currentpt = nextpt;
            } while (k < num[0]);

            vnum[0] = num[0];

            return;
        }

        /* opens files, initialises variables */
        /* prints out the run parameters */
        private void Header()
        {
            max_point = 0;

            string header = "----------HEADER----------\n";

            header += "Dimension: " + D + "\n";
            header += "Volume: " + VOL + "\n";
            header += "Marked node coupling: " + ALPHA + "\n";
            header += "Simplex coupling: " + kappa_d + "\n";
            header += "Bulk coupling: " + BETA + "\n";
            header += "Number of sweeps: " + SWEEPS + "\n";
            header += "Thermalization time: " + THERMALISE + "\n";
            header += "Number of sweeps between KD tuning: " + TUNE_COUPLING + "\n";
            header += "Gap between measurements: " + GAP + "\n";
            header += "Volume fluctuation parameter: " + DV + "\n";
            header += "Random number seed: " + SEED + "\n";

            Debug.Log(header);
        }

        /* opens all the data files and zeroes measure bins */
        private void Init()
        {
            // TODO: is clearing these arrays manually necessary
            for (int i = 0; i < DPLUS; i++)
            {
                try_subsimplex[i] = 0;
                go_subsimplex[i] = 0;
                manifold_subsimplex[i] = 0;
                legal_subsimplex[i] = 0;
            }

            real_simplex = 0f;
            real_node = 0f;
            manifold_check = 0f;
            max_point = 0;
        }

        void InitialConfig()
        {
            int index, tmp;
            int[] dum, dum2;
            dum = new int[DPLUSPLUS];
            dum2 = new int[DPLUSPLUS];

            for (int i = 0; i < DPLUSPLUS; i++)
                dum[i] = i;

            for (int i = 0; i < DPLUSPLUS; i++)
            {
                Combo(dum, dum2, DPLUSPLUS, i);

                simplex_point[pointer_number] = new Simplex(dum2, D, pointer_number);

                simplex_number++;
                pointer_number++;
            }

            /* now set up pointers loop over faces to simplices */

            for (int i = 0; i < DPLUSPLUS; i++)
            {
                for (int j = 0; j < DPLUSPLUS; j++)
                    if (j != i)
                    {
                        index = FindFace(simplex_point[i], dum[j]);
                        simplex_point[i].neighbors[index] = simplex_point[j];
                    }
            }

            string message = "";
            for (int i = 0; i < simplex_number; i++)
            {
                for (int j = 0; j < DPLUS; j++)
                {
                    message += simplex_point[i].vertices[j] + "\t" + simplex_point[i].neighbors[j].label + "\t";
                }
                message += "\n";
            }
            Debug.Log(message);

            node_number = DPLUSPLUS;
            growing_vol = DPLUSPLUS;
            return;
        }

        private void Measure()
        {
            real_simplex += simplex_number;
            real_node += node_number;

            Tidy();

            number_measure++;
        }

        private float DS(int q, int deltaq, float qmean)
        {
            // computes change in local action. Assumes (q-qmean)^2 form
            return deltaq * (deltaq + 2 * q - 2 * qmean);
        }

        private bool Metropolis(int subsimplex, int[] a, Simplex[] addresses)
        {
            float ds0, dsd, dsn, dummy, change;
            bool accept;
            int dum, notnnmarked, imark;
            int[] vnum = new int[1];
            int[] nn = new int[VOL];
            int[] notbound = new int[DPLUSPLUS];
            Simplex p = null;
            int tmp;

            // find marked node
            int done = 0;
            for (int i = 0; i < pointer_number; i++)
            {
                if (simplex_point[i] == null) continue;

                if (done == 1) break;
                for (int j = 0; j < DPLUS; j++)
                {
                    if (simplex_point[i].vertices[j] == 0)
                    {
                        p = simplex_point[i];
                        done = 1;
                        break;
                    }
                }
            }

            // find vertices on boundary
            dum = 0;
            GetNN(p, dum, nn, vnum);

            // flag them
            for (int i = 0; i < DPLUSPLUS; i++)
            {
                notbound[i] = 1;
                for (int j = 0; j < vnum[0]; j++)
                {
                    if (a[i] == nn[j])
                    {
                        notbound[i] = 0;
                        accept = false;
                        return accept;
                    }
                }
            }

            dsn = 0f;
            notnnmarked = 1;
            for (int i = 0; i < DPLUS; i++)
            {
                if (a[i] == 0) notnnmarked = 0;
            }

            ds0 = 0f; // delta N0_boundary=0
            dsd = kappa_d * (2 * subsimplex - D);
            dsd = dsd + (2 * subsimplex - D) * (2 * subsimplex - D + 2 * (simplex_number - vnum[0] - DISKVOL)) / (1f * DV * DV);

            ////////////////////// Node Insertion ////

            if (subsimplex == D)
            {
                /* node insertion */
                for (int i = 0; i < DPLUS; i++)
                {
                    tmp = FindOrder(a[i], addresses[DPLUS]);
                    dsn = dsn + BETA * DS(tmp, D - 1, TC) * notbound[i]; // this only added for bulk points
                }

                //update for the inserted node
                dsn = dsn + BETA * (DPLUS - TC) * (DPLUS - TC) * notnnmarked;
            }
            else if (subsimplex == 0 && a[0] != 0)
            {
                /* node deletion */
                for (int i = 1; i < DPLUS; i++)
                {
                    tmp = FindOrder(a[i], addresses[DPLUS]);
                    dsn = dsn + BETA * DS(tmp, -D + 1, TC) * notbound[i];
                }

                dsn = dsn - BETA * (DPLUS - TC) * (DPLUS - TC) * notnnmarked;
            }

            //////// link flip ///////////////
            else
            {
                for (int i = 0; i < subsimplex; i++)
                {
                    tmp = FindOrder(a[i], addresses[DPLUS]);
                    dsn = dsn + BETA * DS(tmp, 2 * subsimplex - D - 1, TC) * notbound[i];
                }
                for (int i = subsimplex + 1; i < DPLUS; i++)
                {
                    tmp = FindOrder(a[i], addresses[DPLUS]);
                    dsn = dsn + BETA * DS(tmp, 2 * subsimplex - D + 1, TC) * notbound[i];
                }
            }

            if (subsimplex == 0 && a[0] == 0)
            {
                Debug.Log("Trying to delete marked node");
                dsd = 1000000f;
            }

            dsd = dsd + dsn + ds0;
            dsd = Mathf.Exp(dsd);
            dsd = (1f + (2 * subsimplex - D) / (float)(simplex_number)) * dsd;
            dsd = 1f / (1f + dsd);

            dummy = Random.Range(0, 1f) - dsd;

            if (dummy < 0)
                accept = true;
            else
                accept = false;

            if (grow)
                accept = true;

            if (((simplex_number + 2 * subsimplex - D) < MINVOL) || ((simplex_number + 2 * subsimplex - D) > MAXVOL))
                accept = false;

            return accept;
        }

        private void Pop()
        {
            /* if stack is empty do nothing */

            if (stack_head == null)
                return;

            // TODO: just remove temp? we don't return it
            Node temp = stack_head;
            stack_head = stack_head.next;
            stack_count--;

            return;
        }

        /* writes some results to output stream */
        private void PrintOut()
        {
            float dummy = 0f;
            string results = "Results:\n";

            real_simplex = real_simplex / number_measure;
            real_node = real_node / number_measure;

            for (int i = 0; i < DPLUS; i++)
                dummy += legal_subsimplex[i];

            manifold_check /= dummy;
            manifold_check *= VOL;

            float tmp = (real_simplex - VOL) * 2f / (DV * DV);
            results += $"Average number of simplices: {real_simplex}\n";
            results += $"Average number of nodes: {real_node}\n";
            results += $"Average maximum pointer number: {max_point}\n";
            results += $"Average number of simplices in manifold check {manifold_check}\n";
            results += $"Final kappa_d {kappa_d + tmp}\n";

            results += "Subsimplex Moves\n";
            for (int i = 0; i < DPLUS; i++)
            {
                results += $"{i} subsimplices : \n";
                results += $"Number tried {try_subsimplex[i]}\n";
                results += $"Number that are legal {legal_subsimplex[i]}\n";
                results += $"Number that pass manifold test {manifold_subsimplex[i]}\n";
                results += $"Number that pass Metropolis test {go_subsimplex[i]}\n";
            }

            Debug.Log(results);
        }

        /* routine pushes deleted vertex label onto stack - DT.java */
        private void Push(int i)
        {
            Node temp = new Node();
            temp.name = i;
            temp.next = stack_head;

            stack_head = temp;
            stack_count++;
            return;
        }

        /* routine handles reconnection of new simplex pointers */
        private void Reconnect(int[] a, Simplex[] q, int sub)
        {
            int index_face1, index_face2, index_face3, opp;

            /* loop over final state simplices */

            for (int i = 0; i < sub + 1; i++)
            {
                /* 'internal' faces first */

                for (int j = 0; j < sub + 1; j++)
                {
                    if (j != i)
                    {
                        index_face1 = FindFace(q[i], a[j]);
                        q[i].neighbors[index_face1] = q[j];
                    }
                }

                /* now 'external' faces */

                for (int j = sub + 1; j < DPLUSPLUS; j++)
                {
                    index_face1 = FindFace(q[i], a[j]);
                    index_face2 = FindFace(q[j], a[i]);

                    /* have found external simplices involved reconnect outward pointers */

                    q[i].neighbors[index_face1] = q[j].neighbors[index_face2];

                    /* just adjust pointers on external simplices so they point at new ones */

                    opp = q[i].neighbors[index_face1].sum - SumFace(q[i], q[i].vertices[index_face1]);

                    index_face3 = FindFace(q[i].neighbors[index_face1], opp);
                    q[i].neighbors[index_face1].neighbors[index_face3] = q[i];
                }
            }
            return;
        }

        /* selects a simplex at random by accessing an array of pointers */
        /* once has a simplex select subsimplex/move at random */
        private Simplex SelectSimplex(ref int sub)
        {
            Simplex temp;

            do
            {
                int i = Random.Range(0, pointer_number);
                temp = simplex_point[i];
            } while (temp == null);

            /* initially just grow with node insertion moves */

            if (grow)
                sub = D;
            else
                sub = Random.Range(0, DPLUS);

            try_subsimplex[sub]++;

            return temp;
        }

        private void ShiftCoupling()
        {
            float dum = (vol_monitor - b_monitor) / ((float)num_monitor);
            float dum2 = b_monitor / ((float)num_monitor);

            kappa_d = kappa_d + (2f / (DV * DV)) * (dum - DISKVOL) + 2f * ALPHA * (dum2 - MARKEDQ) * FRAC / (2f - FRAC);
            //  ALPHA=ALPHA+2*ALPHA*(dum2-FRAC*VOL/2);

            string message = "";
            message += "New coupling kappa_d is " + kappa_d;
            message += "Average disk volume=" + dum;
            message += "Average boundary length=" + dum2;
            message += "sum of total boundary length=" + b_monitor;
            Debug.Log(message);

            vol_monitor = 0;
            b_monitor = 0;
            num_monitor = 0;

            return;
        }

        /* routine returns sum of vertices around the face conjugate to node i */
        private int SumFace(Simplex p, int i)
        {
            int add = 0;
            for (int j = 0; j < DPLUS; j++)
                if (p.vertices[j] != i) add += p.vertices[j];
            return add;
        }

        /* every sweep clean up pointer array */
        private void Tidy()
        {
            Simplex[] temp = new Simplex[BIGVOL];
            int add = 0;

            /* run down array compressing non NULL extries into new array */
            /* and reassigning simplex labels according to their new index in this */
            /* array. Finally copy back */

            if (pointer_number > max_point) max_point = pointer_number;

            for (int i = 0; i < pointer_number; i++)
            {
                if (simplex_point[i] != null)
                {
                    temp[add] = simplex_point[i];
                    temp[add].label = add;
                    add++;
                }
            }

            for (int i = 0; i < add; i++)
                simplex_point[i] = temp[i];

            pointer_number = add;
            if (pointer_number != simplex_number)
                Debug.Log("oops - pointer number is not equal to simplex_number in tidy()");

            return;
        }

        /* driver for triangulation updates */
        private void TrialChange()
        {
            int subsimplex = D;

            Simplex simp;
            Simplex[] addresses = new Simplex[DPLUSPLUS];
            int[] labels = new int[DPLUSPLUS];
            int[] q = new int[MAXVOL];
            bool legal_move, good_manifold, metro_accept;

            // grab triangle and move type at random
            simp = SelectSimplex(ref subsimplex);

            // check if move is legal i.e coordination of simplex =D+1-subsimplex
            legal_move = GoodSubsimplex(simp, subsimplex, labels, addresses);

            if (!legal_move) return;

            legal_subsimplex[subsimplex]++;

            // make sure move will not create degeneracies
            good_manifold = AllowedMove(simp, subsimplex, labels);

            if (!good_manifold) return;

            manifold_subsimplex[subsimplex]++;

            // check change in action
            metro_accept = Metropolis(subsimplex, labels, addresses);

            if (!metro_accept) return;

            go_subsimplex[subsimplex]++;

            // if accept update triangulation
            DT_Update(labels, addresses, subsimplex);
        }

        /* coordinates addition of new simplices and removal of old ones */
        private void DT_Update(int[] a, Simplex[] q, int sub)
        {
            int[] c = new int[DPLUSPLUS];
            int[] temp = new int[DPLUS];
            /* if subsimplex is node then save its label on the stack */

            if (sub == 0)
            {
                Push(a[0]);
                node_number--;
            }
            if (sub == D)
            {
                Pop();
                node_number++;
            }

            /* loop over new simplices */

            for (int i = 0; i < sub + 1; i++)
            {
                Combo(a, c, sub + 1, i);

                for (int j = sub + 1; j < DPLUSPLUS; j++)
                    c[j - 1] = a[j];

                q[i] = new Simplex(c, D, pointer_number);
                simplex_point[pointer_number] = q[i];

                simplex_number++;
                pointer_number++;
            }

            /* now reconnect pointers appropriately */

            Reconnect(a, q, sub);

            /*  old guys */

            for (int i = sub + 1; i < DPLUSPLUS; i++)
            {
                simplex_point[q[i].label] = null;
                simplex_number--;
            }

            return;
        }
    }
}
