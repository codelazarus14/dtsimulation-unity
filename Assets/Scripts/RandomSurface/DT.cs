using System;
using System.Collections.Generic;
using UnityEngine;
using Random = UnityEngine.Random;

namespace DTSimulation.RandomSurface
{
    public class DT
    {
        // all comments using /**/ are from DT.java

        // TODO: fix style/naming conventions
        // TODO: fix messy/redundant loops and reduce work
        // TODO: remove unused, add cached values, make LOGICS and ints = 0 or 1 become bools (usually "int done"), remove int[1]'s (ref)

        /* some global constants for run */

        private const int DISKVOL = 472;
        private const int MARKEDQ = 282; //FRAC*(DISKVOL+2)/(2.0-FRAC)
        public const int VOL = DISKVOL + MARKEDQ;

        public float Beta = 3f;

        public const int D = 2;
        public const int DPLUS = D + 1;
        public const int DPLUSPLUS = D + 2;

        private const int BIGVOL = 4 * VOL;
        private const int MAXVOL = VOL + 30;

        private const int THERMALISE = 1;

        /* simple global pointers and counters */

        public int node_number = 0;
        public int pointer_number = 0;
        private int simplex_number = 0;

        /* data structures */

        public Simplex[] simplex_point;

        //
        // my fields/deviations from above
        //
        public Vector3[] NodePositions { get; private set; }

        private string configText;

        public DT(string config)
        {
            // TODO: make thermalize a coroutine and set value back up to 500, sweeps to 20000
            configText = config;
        }

        public void Thermalize()
        {
            simplex_point = new Simplex[BIGVOL];
            ReadFile();

            /* thermalise and output info on run */

            Header();

            Debug.Log("Thermalizing lattice");

            for (int therm = 1; therm < THERMALISE; therm++)
            {
                for (int iter = 0; iter < VOL; iter++)
                    TrialChange();
                Tidy();
            }
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
            bool[] mask = new bool[DPLUS];
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
                mask[j] = false;

            for (j = 0; j < n; j++)
                mask[b[j]] = true;

            j = 0;
            for (int i = 0; i < DPLUS; i++)
            {
                if (!mask[i])
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

                a[DPLUS] = node_number;

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

            // load config
            string[] configLines = configText.Split(new char[0], StringSplitOptions.RemoveEmptyEntries);
            // using a string stack to imitate use of the Scanner in the java version
            Stack<string> configStack = new Stack<string>();
            for (int i = configLines.Length - 1; i >= 0; i--)
                configStack.Push(configLines[i]);

            // read config values
            try
            {
                int simplexCount = int.Parse(configStack.Pop());
                int[,] neighborBuffer = new int[simplexCount, DPLUS];
                int[] vertexBuffer = new int[DPLUS];

                node_number = int.Parse(configStack.Pop());
                int stackCount = int.Parse(configStack.Pop());
                // TODO: why do these exist + debug cleanup
                float.Parse(configStack.Pop());
                float.Parse(configStack.Pop());
                float.Parse(configStack.Pop());

                Debug.Log($"{simplexCount}, {node_number}, {stackCount}");

                simplex_number = 0;
                pointer_number = 0;

                Debug.Log("Reading in existing configuration");

                Debug.Log($"Configuration has volume: {simplexCount}");
                Debug.Log($"Node number: {node_number}");
                Debug.Log($"curvaturesq coupling: {Beta}");

                for (int i = 0; i < stackCount; i++)
                {
                    // ignore stack labels for now
                    int.Parse(configStack.Pop());
                    //// add used vertex labels to stack
                    //dummy = int.Parse(configStack.Pop());
                    //Push(dummy);
                }

                // set up simplices
                for (int i = 0; i < simplexCount; i++)
                {
                    for (int j = 0; j < DPLUS; j++)
                    {
                        vertexBuffer[j] = int.Parse(configStack.Pop());
                        neighborBuffer[i, j] = int.Parse(configStack.Pop());
                    }

                    simplex_point[i] = new Simplex(vertexBuffer, DPLUS, pointer_number);
                    simplex_number++;
                    pointer_number++;
                }

                // set neighbors
                for (int i = 0; i < simplexCount; i++)
                {
                    for (int j = 0; j < DPLUS; j++)
                    {
                        simplex_point[i].neighbors[j] = simplex_point[neighborBuffer[i, j]];
                    }
                }

                // give nodes random positions
                NodePositions = new Vector3[node_number];
                for (int i = 0; i < node_number; i++)
                {
                    NodePositions[i] = new Vector3(Random.Range(-1f, 1f), Random.Range(-1f, 1f), Random.Range(-1f, 1f));
                }

                // print out results
                string message = "";
                for (int i = 0; i < simplexCount; i++)
                {
                    for (int j = 0; j < DPLUS; j++)
                    {
                        message += simplex_point[i].vertices[j] + "\t" + simplex_point[i].neighbors[j].label + "\t";
                    }
                    message += "\n";
                }
                Debug.Log(message);

                Debug.Log("Have read data successfully");
            }
            catch (FormatException e)
            {
                Debug.LogError($"Error parsing config file! {e.Message}");
            }
        }

        // TODO: remove or squash with GetNN and OrderedNeighbors
        // GetNN for nodes
        public (int[], int) NearestNeighbors(int nodeLabel)
        {
            Simplex simp = null;

            for (int i = 0; i < pointer_number && simp == null; i++)
            {
                if (simplex_point[i] == null) continue;
                for (int j = 0; j < DPLUS; j++)
                {
                    if (simplex_point[i].vertices[j] == nodeLabel)
                    {
                        simp = simplex_point[i];
                        break;
                    }
                }
            }

            if (simp == null)
                Debug.LogError($"couldn't find neighbors for node {nodeLabel}");

            int[] neighbors = new int[VOL], nCount = new int[1];
            NearestNeighbors(simp, nodeLabel, neighbors, nCount);
            return (neighbors, nCount[0]);
        }

        // hardwired for D=2 right now ..
        private void NearestNeighbors(Simplex simp, int node, int[] neighbors, int[] nCount)
        {
            Simplex[] simplices = new Simplex[BIGVOL];
            int[] dummy = new int[DPLUS];
            int[] v1 = new int[VOL];
            int[] v2 = new int[VOL];
            bool[] seen = new bool[VOL];
            int k;

            dummy[0] = node;
            FindSimplices(simp, dummy, 1, simplices, nCount);

            k = 0;
            for (int i = 0; i < nCount[0]; i++)
            {
                int index = 0;

                for (int j = 0; j < DPLUS; j++)
                {
                    // find center for neighbor search
                    if (simplices[i].vertices[j] == node)
                    {
                        index = j;
                        break;
                    }
                }

                // find other two points on simplex
                v1[k] = simplices[i].vertices[(index + 1) % DPLUS];
                v2[k] = simplices[i].vertices[(index + 2) % DPLUS];
                k++;
            }

            // mark center node as seen
            neighbors[0] = v1[0];
            seen[0] = true;
            k = 1;

            // look for v1[0] in rest of v1/v2 arrays
            int currNode = v1[0];
            do
            {
                int nextNode = 0;
                for (int i = 0; i < nCount[0]; i++)
                {
                    if (!seen[i])
                    {
                        if (v1[i] == currNode)
                        {
                            neighbors[k] = v2[i];
                            nextNode = v2[i];
                            seen[i] = true;
                            k++;
                        }
                        if (v2[i] == currNode)
                        {
                            neighbors[k] = v1[i];
                            nextNode = v1[i];
                            seen[i] = true;
                            k++;
                        }
                    }
                }
                currNode = nextNode;
            } while (k < nCount[0]);

            return;
        }

        public void RelabelNodes()
        {
            const int VERYBIG = 100000;
            int new2, old;

            int[] num = new int[1];

            bool[] not_seen = new bool[VOL];
            Simplex[] list = new Simplex[VOL];
            int[] dum = new int[DPLUS];

            for (int i = 0; i < VOL; i++)
                not_seen[i] = true;

            new2 = VERYBIG + 1;

            for (int i = 0; i < simplex_number; i++)
            {
                for (int j = 0; j < DPLUS; j++)
                {
                    old = simplex_point[i].vertices[j];

                    if (old == 0 || old > VERYBIG) continue;
                    if (not_seen[old])
                    {
                        dum[0] = old;
                        int dummy = 1;
                        FindSimplices(simplex_point[i], dum, dummy, list, num);
                        for (int k = 0; k < num[0]; k++)
                        {
                            for (int l = 0; l < DPLUS; l++)
                            {
                                if (list[k].vertices[l] == old)
                                    list[k].vertices[l] = new2;
                            }
                        }
                        not_seen[old] = false;
                        new2++;
                    }
                }
            }

            // reset labels
            int a, add;
            for (int i = 0; i < simplex_number; i++)
            {
                add = 0;
                for (int j = 0; j < DPLUS; j++)
                {
                    a = simplex_point[i].vertices[j];
                    a = a % VERYBIG;
                    simplex_point[i].vertices[j] = a;
                    add += a;
                }
                simplex_point[i].sum = add;
            }
        }

        /* opens files, initialises variables */
        /* prints out the run parameters */
        private void Header()
        {
            string header = "----------HEADER----------\n";

            header += "Dimension: " + D + "\n";
            header += "Volume: " + VOL + "\n";
            header += "Bulk coupling: " + Beta + "\n";
            header += "Thermalization time: " + THERMALISE + "\n";

            Debug.Log(header);
        }

        private float DS(Vector3[] vPosns)
        {
            float kl = (vPosns[3] - vPosns[2]).sqrMagnitude;
            float ij = (vPosns[1] - vPosns[0]).sqrMagnitude;

            return kl - ij;
        }

        private float DS_Curvature(Vector3[] posns, Vector3[] vPosns)
        {

            // calculate normals (counterclockwise)
            // for labels [i, j, k, l] = indices [0, 3]
            //  v2 --- k ---- v1
            //   \   / | \   /
            //    \ /  |  \ / 
            //     i --|-- j
            //    / \  |  / \
            //   /   \ | /   \
            //  v4 --- l --- v3

            // before link flip
            // n1 = ij X ik
            // n2 = il X ij
            Vector3 innerBefore1 = Vector3.Cross(posns[1] - posns[0], posns[2] - posns[0]).normalized;
            Vector3 innerBefore2 = Vector3.Cross(posns[3] - posns[0], posns[1] - posns[0]).normalized;
            // after link flip
            // n1 = lk X li
            // n2 = lj X lk
            Vector3 innerAfter1 = Vector3.Cross(posns[2] - posns[3], posns[0] - posns[3]).normalized;
            Vector3 innerAfter2 = Vector3.Cross(posns[1] - posns[3], posns[2] - posns[3]).normalized;

            // triangle neighbors of vertices along the outer edges (i-k-j-l-i around)
            // n1 = jv1 X jk
            // n2 = ik X iv2
            // n3 = jl X jv3
            // n4 = iv4 X il
            Vector3 outer1 = Vector3.Cross(vPosns[0] - posns[1], posns[2] - posns[1]).normalized;
            Vector3 outer2 = Vector3.Cross(posns[2] - posns[0], vPosns[1] - posns[0]).normalized;
            Vector3 outer3 = Vector3.Cross(posns[3] - posns[1], vPosns[2] - posns[1]).normalized;
            Vector3 outer4 = Vector3.Cross(vPosns[3] - posns[0], posns[3] - posns[0]).normalized;

            float beforeDots =
                Vector3.Dot(innerBefore1, innerBefore2) +
                Vector3.Dot(outer1, innerBefore1) +
                Vector3.Dot(outer2, innerBefore1) +
                Vector3.Dot(outer3, innerBefore2) +
                Vector3.Dot(outer4, innerBefore2);
            float afterDots =
                Vector3.Dot(innerAfter1, innerAfter2) +
                Vector3.Dot(outer1, innerAfter1) +
                Vector3.Dot(outer2, innerAfter1) +
                Vector3.Dot(outer3, innerAfter2) +
                Vector3.Dot(outer4, innerAfter2);

            return -Beta * (afterDots - beforeDots);
        }

        private int[] GetOuterSimplexNeighbors(int[] labels, Simplex[] addresses)
        {
            // hardcoded for 2D 
            // two triangles (2,3) and list of the four vertices (i/j shared)
            // find simplex from findface(2/3, i/j) - 4 total

            // link[1, 2] --> v1
            // link[0, 2] --> v2
            // link[1, 3] --> v3
            // link[0, 3] --> v4
            int[] result = new int[4];

            for (int i = 0; i < 4; i++)
            {
                Simplex p = addresses[i / 2 + 2];   // addresses[2/3]
                int center = labels[i % 2];         // labels[i/j]
                Simplex n = p.neighbors[FindFace(p, center)];
                // get label of the vertex not included in the sum (opposite the link)
                result[i] = n.sum - SumFace(p, center);
            }

            return result;
        }

        private bool Metropolis(int subsimplex, int[] a, Simplex[] addresses)
        {
            // TODO: un-fix the volume

            //////// link flip ///////////////

            // get node positions
            Vector3[] posns = new Vector3[a.Length];

            int[] outerNeighbors = GetOuterSimplexNeighbors(a, addresses);
            Vector3[] oPosns = new Vector3[outerNeighbors.Length];

            for (int i = 0; i < a.Length; i++)
            {
                posns[i] = NodePositions[a[i]];
                oPosns[i] = NodePositions[outerNeighbors[i]];
            }

            // deltaE = link flip change + change in curvature
            float dsd = DS(posns) + DS_Curvature(posns, oPosns);

            return MetropolisTest(dsd);
        }

        private bool MetropolisTest(float deltaE)
        {
            deltaE = Mathf.Exp(-deltaE);
            float dummy = Random.Range(0, 1f) - deltaE;

            return dummy < 0;
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

            //sub = Random.Range(0, DPLUS);
            // force only link flips for now
            sub = 1;

            return temp;
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
        public void Tidy()
        {
            Simplex[] temp = new Simplex[BIGVOL];
            int add = 0;

            /* run down array compressing non NULL extries into new array */
            /* and reassigning simplex labels according to their new index in this */
            /* array. Finally copy back */

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
        public void TrialChange()
        {
            int subsimplex = D;

            Simplex simp;
            Simplex[] addresses = new Simplex[DPLUSPLUS];
            int[] labels = new int[DPLUSPLUS];
            bool legal_move, good_manifold, metro_accept;

            // grab triangle and move type at random
            simp = SelectSimplex(ref subsimplex);

            // check if move is legal i.e coordination of simplex =D+1-subsimplex
            legal_move = GoodSubsimplex(simp, subsimplex, labels, addresses);

            if (!legal_move) return;

            // make sure move will not create degeneracies
            good_manifold = AllowedMove(simp, subsimplex, labels);

            if (!good_manifold) return;

            // check change in action
            metro_accept = Metropolis(subsimplex, labels, addresses);

            if (!metro_accept) return;

            // if accept update triangulation
            DT_Update(labels, addresses, subsimplex);
        }

        /* coordinates addition of new simplices and removal of old ones */
        private void DT_Update(int[] a, Simplex[] q, int sub)
        {
            int[] c = new int[DPLUSPLUS];
            int[] temp = new int[DPLUS];

            /* loop over new simplices */

            for (int i = 0; i < sub + 1; i++)
            {
                Combo(a, c, sub + 1, i);

                for (int j = sub + 1; j < DPLUSPLUS; j++)
                    c[j - 1] = a[j];

                q[i] = new Simplex(c, DPLUS, pointer_number);
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

        public void WobbleVertex()
        {
            // size of wobble
            const float magnitude = 0.1f;

            Simplex randSimplex;
            do
            {
                int i = Random.Range(0, pointer_number);
                randSimplex = simplex_point[i];
            } while (randSimplex == null);

            int randNode = randSimplex.vertices[Random.Range(0, DPLUS)];

            Vector3 pos = NodePositions[randNode];

            // create new pos with random wobble
            Vector3 wobble = new Vector3(Random.Range(-1f, 1f), Random.Range(-1f, 1f), Random.Range(-1f, 1f)).normalized;
            Vector3 newPos = pos + wobble * magnitude;

            // find neighbor verts + outer verts of simplices bordering inner group of triangles
            (Simplex[] nSimplices, int[] nNodes) = GetOrderedNeighbors(randSimplex, randNode);

            int[] outerNeighbors = new int[nNodes.Length];

            // TODO: duplicated from GetOuterSimplexNeighbors
            for (int i = 0; i < outerNeighbors.Length; i++)
            {
                Simplex p = nSimplices[i];
                Simplex n = p.neighbors[FindFace(p, randNode)];
                // get label of the vertex not included in the sum (opposite the link)
                outerNeighbors[i] = n.sum - SumFace(p, randNode);
            }

            // compute difference in energy
            float deltaE = DS_WobbleLinks(pos, newPos, nNodes) + DS_WobbleCurvature(pos, newPos, nNodes, outerNeighbors);

            // update position with wobble if accepted
            if (MetropolisTest(deltaE))
            {
                NodePositions[randNode] = newPos;
            }
        }

        private float DS_WobbleLinks(Vector3 before, Vector3 after, int[] neighbors)
        {
            float beforeTotal = 0f;
            float afterTotal = 0f;

            for (int i = 0; i < neighbors.Length; i++)
            {
                // get total distances from each node to center, before and after
                beforeTotal += (NodePositions[neighbors[i]] - before).sqrMagnitude;
                afterTotal += (NodePositions[neighbors[i]] - after).sqrMagnitude;
            }

            return afterTotal - beforeTotal;
        }

        private float DS_WobbleCurvature(Vector3 before, Vector3 after, int[] neighbors, int[] outerNeighbors)
        {
            // compute difference in normals
            return -Beta * (GetCurvature(after, neighbors, outerNeighbors) - GetCurvature(before, neighbors, outerNeighbors));
        }

        private float GetCurvature(Vector3 centerPos, int[] neighbors, int[] outerNeighbors)
        {
            // collect all the normals of the surrounding simplices
            Vector3[] normals = new Vector3[neighbors.Length];
            // normals of the triangles sharing boundary links [p1-p2], [p2-p3] etc.
            Vector3[] outerNormals = new Vector3[outerNeighbors.Length];

            for (int i = 0; i < neighbors.Length - 1; i++)
            {
                int p1 = neighbors[i];
                int p2 = neighbors[i + 1];
                int v = outerNeighbors[i];

                normals[i] = Vector3.Cross(NodePositions[p1] - centerPos, NodePositions[p2] - centerPos).normalized;
                outerNormals[i] = Vector3.Cross(NodePositions[v] - NodePositions[p1], NodePositions[p2] - NodePositions[p1]).normalized;
            }
            normals[^1] = Vector3.Cross(NodePositions[neighbors[^1]] - centerPos, NodePositions[neighbors[0]] - centerPos).normalized;
            outerNormals[^1] = Vector3.Cross(NodePositions[outerNeighbors[^1]] - NodePositions[neighbors[^1]], NodePositions[neighbors[0]] - NodePositions[neighbors[^1]]).normalized;

            // dot each of the pairs (prev-next inner and inner-outer)
            float curv = 0f;
            for (int i = 0; i < neighbors.Length - 1; i++)
            {
                curv += Vector3.Dot(normals[i], normals[i + 1]);
                curv += Vector3.Dot(normals[i], outerNormals[i]);
            }
            curv += Vector3.Dot(normals[^1], normals[0]);
            curv += Vector3.Dot(normals[^1], outerNormals[^1]);

            return curv;
        }

        private (Simplex[] simplices, int[] nodes) GetOrderedNeighbors(Simplex s, int center)
        {
            List<Simplex> simplices = new List<Simplex>();
            List<int> nodes = new List<int>();

            // simplex currently being evaluated
            Simplex currSimplex = s;
            // used with FindSimplices to pick one of D simplices sharing a link
            Simplex[] nearSimplices = new Simplex[D];
            int[] nCount = new int[1];
            int[] link = new int[D];

            do
            {
                // set up link
                link[0] = center;
                foreach (int vertex in currSimplex.vertices)
                {
                    if (!(vertex == center || (nodes.Count > 0 && vertex == nodes[^1])))
                    {
                        link[1] = vertex;
                        break;
                    }
                }

                // get two triangles that share the link
                FindSimplices(currSimplex, link, 2, nearSimplices, nCount);

                // find the (unique) new simplex
                foreach (Simplex simp in nearSimplices)
                {
                    if (simp != currSimplex)
                    {
                        currSimplex = simp;
                        break;
                    }
                }
                // add new neighbor to lists
                nodes.Add(link[1]);
                simplices.Add(currSimplex);
            } while (currSimplex != s);

            return (simplices.ToArray(), nodes.ToArray());
        }
    }
}
