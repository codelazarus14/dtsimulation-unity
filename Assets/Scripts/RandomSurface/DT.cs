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

        public DT(DTConfig cfg)
        {
            // TODO: make thermalize a coroutine and set value back up to 500, sweeps to 20000
            LoadConfig(cfg);
        }

        public void Thermalize()
        {
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
            Simplex[] examined = new Simplex[VOL];
            Simplex[] remaining = new Simplex[VOL];
            Simplex[] newNeighbors = new Simplex[VOL];
            Simplex[] common = new Simplex[DPLUS];
            int[] commonVerts = new int[DPLUS];
            int numRemaining = 1;
            int numExamined = 1;
            bool good = true;

            /* need to teminate when either 1. see opposing vertex a[DPLUS] in some */
            /* simplex or 2. have examined all simplices containing the new common */
            /* vertices derived from original simplex */
            /* so need to flag simplices that have been examined */

            if (sub == D || sub == 0) return true;

            for (int i = sub + 1; i < DPLUS; i++)
                commonVerts[i - (sub + 1)] = a[i];

            examined[0] = p;
            remaining[0] = p;
            p.seen = true;

            /* loop while new neighbour simplices to examine */

            while (numRemaining > 0)
            {
                int newNeighborCount = 0;

                for (int i = 0; i < numRemaining; i++)
                {
                    /* examine to see if contains 'opposing vertex */

                    for (int j = 0; j < DPLUS && good; j++)
                    {
                        if (remaining[i].vertices[j] == a[DPLUS])
                            good = false;
                    }

                    if (!good) break;

                    /* find simplices which also in common with this subsimplex */
                    /* and neighbour to simplex array1[i] */

                    CommonSimplices(remaining[i], commonVerts, D - sub, ref common);

                    /* check to see whether any of these have been seen before */
                    /* if not then add to array2 and flag seen. Also make note in examine */

                    for (int j = 0; j < sub + 1; j++)
                    {
                        if (!common[j].seen)
                        {
                            newNeighbors[newNeighborCount++] = common[j];
                            common[j].seen = true;

                            examined[numExamined++] = common[j];
                        }
                    }
                }

                if (!good) break;

                for (int i = 0; i < newNeighborCount; i++)
                    remaining[i] = newNeighbors[i];

                numRemaining = newNeighborCount;
            }

            /* depending on value of good we have either found the opposing vertex or */
            /* examined all the associated simplices */

            /* first set all local flags back to zero */

            for (int i = 0; i < numExamined; i++)
                examined[i].seen = false;

            return good;
        }

        /* routine takes a named subsimplex a[0]...a[n-1] and searches for it */
        /* in a simplex pointed at by p. */
        /* returns a vector containing the local indices of pointers to */
        /* neighbouring simplices */
        /* which share a face also encompassing this subsimplex */
        private void CommonSimplices(Simplex p, int[] a, int n, ref Simplex[] face)
        {
            const int NOTFOUND = -1;
            int[] commonIndices = new int[n];
            bool[] mask = new bool[DPLUS];

            for (int i = 0; i < n; i++)
                commonIndices[i] = NOTFOUND;

            /* find positions/local indices of subsimplex in this simplex */

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < DPLUS && commonIndices[i] == NOTFOUND; j++)
                {
                    if (p.vertices[j] == a[i])
                        commonIndices[i] = j;
                }

                if (commonIndices[i] == NOTFOUND)
                {
                    Debug.LogError("Error in CommonSimplex");
                    return;
                }
                mask[commonIndices[i]] = true;
            }

            int commonIdx = 0;
            for (int i = 0; i < DPLUS; i++)
            {
                if (!mask[i])
                    face[commonIdx++] = p.neighbors[i];
            }
        }

        /* routine takes pointer to simplex and a vertex and returns local */
        /* index to conjugate face */
        private int FindFace(Simplex p, int a)
        {
            for (int i = 0; i < DPLUS; i++)
                if (p.vertices[i] == a) return i;

            /* return dummy if get to here error */
            Debug.LogError("Error in FindFace");
            return VOL;
        }

        /* finds addresses of all simplices which share a given subsimplex */
        private void FindSimplices(Simplex p, int[] a, int sub, ref Simplex[] simplices, out int numSimplices)
        {
            // TODO: a lot of this code seems to be copied from AllowedMove
            Simplex[] remaining = new Simplex[VOL];
            Simplex[] newNeighbors = new Simplex[VOL];
            Simplex[] common = new Simplex[DPLUS];
            int numRemaining = 1;
            int numSeen = 1;

            remaining[0] = p;
            simplices[0] = p;
            p.seen = true;

            while (numRemaining > 0)
            {
                int newNeighborCount = 0;

                for (int i = 0; i < numRemaining; i++)
                {
                    CommonSimplices(remaining[i], a, sub, ref common);

                    for (int j = 0; j < DPLUS - sub; j++)
                    {
                        if (!common[j].seen)
                        {
                            simplices[numSeen++] = common[j];
                            newNeighbors[newNeighborCount++] = common[j];
                            common[j].seen = true;
                        }
                    }
                }

                for (int i = 0; i < newNeighborCount; i++)
                    remaining[i] = newNeighbors[i];

                numRemaining = newNeighborCount;
            }

            for (int i = 0; i < numSeen; i++)
                simplices[i].seen = false;

            numSimplices = numSeen;
        }

        /* takes simplex and order of subsimplex returns logical flag if legal move */
        /* legal move is equivalent to being exactly d+1-i associated simplices. */
        /* these subsimplex vertices occupy first i+1 entries in a[] */
        /* other d-i 'external' vertices of original simplex occupy end of a[] */
        /* also returns pointers to these d+1-i simplices */
        /* and opposing vertex is placed at end of a[] */
        private bool GoodSubsimplex(Simplex p, int sub, ref int[] a, ref Simplex[] addresses)
        {
            int opposing;
            int[] aind = new int[DPLUS];

            /* generate subsimplex at random placing its indices in aind */

            int[] randIndices = new int[DPLUS];

            for (int i = 0; i < DPLUS; i++)
                randIndices[i] = i;

            // shuffle to get random indices
            int currIdx = DPLUS;
            while (currIdx > 1)
            {
                int k = Random.Range(0, currIdx--);
                (randIndices[k], randIndices[currIdx]) = (randIndices[currIdx], randIndices[k]);
            }

            // copy random indices into aind
            for (int i = 0; i < DPLUS; i++)
                aind[i] = randIndices[i];

            /* now loop over all possible faces constucted to include this subsimplex */
            /* by selecting d-i-1 out of the d-i remaining indices */

            for (int i = 0; i < DPLUS; i++)
                a[i] = p.vertices[aind[i]];

            addresses[sub + 1] = p.neighbors[aind[sub + 1]];
            opposing = addresses[sub + 1].sum - SumFace(p, a[sub + 1]);

            /* protect following loop if sub=D-1 */

            if (sub < (D - 1))
            {
                for (int i = sub + 2; i < DPLUS; i++)
                {
                    addresses[i] = p.neighbors[aind[i]];
                    if (addresses[i].sum - SumFace(p, a[i]) != opposing)
                        return false;
                }
            }

            a[DPLUS] = opposing;
            addresses[DPLUS] = p;
            return true;
        }

        private void LoadConfig(DTConfig config)
        {
            simplex_point = new Simplex[BIGVOL];
            simplex_number = config.simplexCount;
            node_number = config.nodeCount;

            int[,] neighborBuffer = new int[simplex_number, DPLUS];
            int[] vertexBuffer = new int[DPLUS];

            // set up simplices
            for (int i = 0; i < simplex_number; i++)
            {
                for (int j = 0; j < DPLUS; j++)
                {
                    vertexBuffer[j] = config.vertices[i * DPLUS + j];
                    neighborBuffer[i, j] = config.neighbors[i * DPLUS + j];
                }

                simplex_point[i] = new Simplex(vertexBuffer, DPLUS, pointer_number++);
            }

            // set neighbors
            for (int i = 0; i < simplex_number; i++)
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
            for (int i = 0; i < simplex_number; i++)
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

        // GetNN for nodes
        public int[] NearestNeighbors(int nodeLabel)
        {
            Simplex p = null;

            // find a simplex which contains the node
            for (int i = 0; i < pointer_number && p == null; i++)
            {
                if (simplex_point[i] == null) continue;
                for (int j = 0; j < DPLUS; j++)
                {
                    if (simplex_point[i].vertices[j] == nodeLabel)
                    {
                        p = simplex_point[i];
                        break;
                    }
                }
            }

            if (p == null)
                Debug.LogError($"couldn't find neighbors for node {nodeLabel}");

            (_, int[] neighbors) = GetOrderedNeighbors(p, nodeLabel);

            return neighbors;
        }

        public void RelabelNodes()
        {
            // number larger than any node label
            const int VERYBIG = 100000;
            int newLabel = VERYBIG + 1, oldLabel;

            Simplex[] nearSimplices = new Simplex[VOL];
            bool[] seen = new bool[VOL];
            int[] dum = new int[DPLUS];

            for (int i = 0; i < simplex_number; i++)
            {
                for (int j = 0; j < DPLUS; j++)
                {
                    oldLabel = simplex_point[i].vertices[j];

                    // skip marked node and any already-updated labels
                    if (oldLabel == 0 || oldLabel > VERYBIG) continue;

                    if (!seen[oldLabel])
                    {
                        // loop over all simplices that share this node
                        dum[0] = oldLabel;
                        FindSimplices(simplex_point[i], dum, 1, ref nearSimplices, out int nCount);

                        for (int k = 0; k < nCount; k++)
                        {
                            for (int l = 0; l < DPLUS; l++)
                            {
                                // update vertices with new label
                                if (nearSimplices[k].vertices[l] == oldLabel)
                                    nearSimplices[k].vertices[l] = newLabel;
                            }
                        }
                        seen[oldLabel] = true;
                        newLabel++;
                    }
                }
            }

            // reset labels, resulting nodes labeled 0-node_number
            for (int i = 0; i < simplex_number; i++)
            {
                int finalSum = 0;
                for (int j = 0; j < DPLUS; j++)
                {
                    int finalLabel = simplex_point[i].vertices[j];
                    finalLabel = finalLabel % VERYBIG;
                    simplex_point[i].vertices[j] = finalLabel;
                    finalSum += finalLabel;
                }
                // update simplex with sum of new labels
                simplex_point[i].sum = finalSum;
            }
        }

        /* opens files, initialises variables */
        /* prints out the run parameters */
        private void Header()
        {
            string header = "----------HEADER----------\n";

            header += $"Dimension: {D}\n";
            header += $"Volume: {VOL}\n";
            header += $"Bulk coupling: {Beta}\n";
            header += $"Thermalization time: {THERMALISE}\n";

            Debug.Log(header);
        }

        private float DS(Vector3[] vPosns)
        {
            // return change in energy between the two links (before and after flip)
            float kl = (vPosns[3] - vPosns[2]).sqrMagnitude;
            float ij = (vPosns[1] - vPosns[0]).sqrMagnitude;

            return kl - ij;
        }

        private float DS_Curvature(Vector3[] posns, Vector3[] vPosns)
        {
            // calculate normals (counterclockwise)
            // see GetOuterSimplexNeighbors for diagram

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
            // labels is a list of the four vertices (i/j shared), addresses has two triangles (a[2],a[3])

            // for labels [i, j, k, l] = indices [0, 3]
            //  v2 --- k ---- v1
            //   \   / | \   /
            //    \ /  |  \ / 
            //     i --|-- j
            //    / \  |  / \
            //   /   \ | /   \
            //  v4 --- l --- v3

            // link[1, 2] --> v1
            // link[0, 2] --> v2
            // link[1, 3] --> v3
            // link[0, 3] --> v4
            int[] result = new int[4];

            for (int i = 0; i < 4; i++)
            {
                Simplex p = addresses[i / 2 + 2];   // addresses[2/3]
                int origin = labels[i % 2];         // labels[i/j]
                // get neighbor of origin index within each simplex - outer simplices
                Simplex n = p.neighbors[FindFace(p, origin)];
                // get label of the outer simplex vertex not included in the sum (opposite the link)
                result[i] = n.sum - SumFace(p, origin);
            }

            return result;
        }

        private bool Metropolis(int _, int[] a, Simplex[] addresses)
        {
            //////// link flip ///////////////     only when sub = 1

            // get node positions
            Vector3[] posns = new Vector3[a.Length];

            int[] outerNeighbors = GetOuterSimplexNeighbors(a, addresses);
            Vector3[] outerPosns = new Vector3[outerNeighbors.Length];

            for (int i = 0; i < a.Length; i++)
            {
                posns[i] = NodePositions[a[i]];
                outerPosns[i] = NodePositions[outerNeighbors[i]];
            }

            // deltaE = link flip change + change in curvature
            float dsd = DS(posns) + DS_Curvature(posns, outerPosns);

            return MetropolisTest(dsd);
        }

        private bool MetropolisTest(float deltaE)
        {
            deltaE = Mathf.Exp(-deltaE);
            float dummy = Random.Range(0, 1f) - deltaE;

            return dummy < 0;
        }

        /* routine handles reconnection of new simplex pointers */
        private void Reconnect(int[] a, ref Simplex[] simplices, int sub)
        {
            /* loop over final state simplices */

            for (int i = 0; i < sub + 1; i++)
            {
                /* 'internal' faces first */

                for (int j = 0; j < sub + 1; j++)
                {
                    if (j != i)
                    {
                        int face1Index = FindFace(simplices[i], a[j]);
                        simplices[i].neighbors[face1Index] = simplices[j];
                    }
                }

                /* now 'external' faces */

                for (int j = sub + 1; j < DPLUSPLUS; j++)
                {
                    int face1Index = FindFace(simplices[i], a[j]);
                    int face2Index = FindFace(simplices[j], a[i]);

                    /* have found external simplices involved reconnect outward pointers */

                    simplices[i].neighbors[face1Index] = simplices[j].neighbors[face2Index];

                    /* just adjust pointers on external simplices so they point at new ones */

                    int opposing = simplices[i].neighbors[face1Index].sum - SumFace(simplices[i], simplices[i].vertices[face1Index]);

                    int face3Index = FindFace(simplices[i].neighbors[face1Index], opposing);
                    simplices[i].neighbors[face1Index].neighbors[face3Index] = simplices[i];
                }
            }
        }

        /* selects a simplex at random by accessing an array of pointers */
        /* once has a simplex select subsimplex/move at random */
        private Simplex SelectSimplex(out int sub)
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
            int vertexSum = 0;
            for (int j = 0; j < DPLUS; j++)
                if (p.vertices[j] != i) vertexSum += p.vertices[j];
            return vertexSum;
        }

        /* every sweep clean up pointer array */
        public void Tidy()
        {
            Simplex[] temp = new Simplex[BIGVOL];
            int newSimplexCount = 0;

            /* run down array compressing non NULL extries into new array */
            /* and reassigning simplex labels according to their new index in this */
            /* array. Finally copy back */

            for (int i = 0; i < pointer_number; i++)
            {
                if (simplex_point[i] != null)
                {
                    temp[newSimplexCount] = simplex_point[i];
                    temp[newSimplexCount].label = newSimplexCount;
                    newSimplexCount++;
                }
            }

            for (int i = 0; i < newSimplexCount; i++)
                simplex_point[i] = temp[i];

            pointer_number = newSimplexCount;
            if (pointer_number != simplex_number)
                Debug.Log("oops - pointer number is not equal to simplex_number in tidy()");
        }

        /* driver for triangulation updates */
        public void TrialChange()
        {
            Simplex[] addresses = new Simplex[DPLUSPLUS];
            int[] labels = new int[DPLUSPLUS];
            bool isLegal, isGoodManifold, metroAccepts;

            // grab triangle and move type at random
            Simplex p = SelectSimplex(out int subsimplex);

            // check if move is legal i.e coordination of simplex =D+1-subsimplex
            isLegal = GoodSubsimplex(p, subsimplex, ref labels, ref addresses);

            if (!isLegal) return;

            // make sure move will not create degeneracies
            isGoodManifold = AllowedMove(p, subsimplex, labels);

            if (!isGoodManifold) return;

            // check change in action
            metroAccepts = Metropolis(subsimplex, labels, addresses);

            if (!metroAccepts) return;

            // if accept update triangulation
            Update(labels, ref addresses, subsimplex);
        }

        /* routine takes n element vector a[] and returns n-1 element in b[] */
        /* thus n combinations possible selected according to count */
        /* count is the index 'left out' in getting n-1 from n */
        /* vertex left out of final vector b[] is returned in b[n-1] */
        private void Combo(int[] a, ref int[] b, int n, int leaveOut)
        {
            int currIdx = 0;
            for (int i = 0; i < n; i++)
            {
                if (leaveOut != i)
                {
                    b[currIdx++] = a[i];
                }
                else b[n - 1] = a[i];
            }
        }

        /* coordinates addition of new simplices and removal of old ones */
        private void Update(int[] a, ref Simplex[] newSimplices, int sub)
        {
            int[] newVerts = new int[DPLUSPLUS];

            /* loop over new simplices */

            for (int i = 0; i < sub + 1; i++)
            {
                Combo(a, ref newVerts, sub + 1, i);

                for (int j = sub + 1; j < DPLUSPLUS; j++)
                    newVerts[j - 1] = a[j];

                newSimplices[i] = new Simplex(newVerts, DPLUS, pointer_number);
                simplex_point[pointer_number++] = newSimplices[i];

                simplex_number++;
            }

            /* now reconnect pointers appropriately */

            Reconnect(a, ref newSimplices, sub);

            /*  old guys */

            for (int i = sub + 1; i < DPLUSPLUS; i++)
            {
                simplex_point[newSimplices[i].label] = null;
                simplex_number--;
            }
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

        private (Simplex[] simplices, int[] nodes) GetOrderedNeighbors(Simplex p, int center)
        {
            List<Simplex> simplices = new List<Simplex>();
            List<int> nodes = new List<int>();

            // simplex currently being evaluated
            Simplex currSimplex = p;
            // used with FindSimplices to pick one of D simplices sharing a link
            Simplex[] nearSimplices = new Simplex[D];
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
                FindSimplices(currSimplex, link, 2, ref nearSimplices, out int nCount);

                // find the (unique) new simplex
                for (int i = 0; i < nCount; i++)
                {
                    if (nearSimplices[i] != currSimplex)
                    {
                        currSimplex = nearSimplices[i];
                        break;
                    }
                }
                // add new neighbor to lists
                nodes.Add(link[1]);
                simplices.Add(currSimplex);
            } while (currSimplex != p);

            return (simplices.ToArray(), nodes.ToArray());
        }
    }
}
