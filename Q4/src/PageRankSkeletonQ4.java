 /*
  * Use command-line flag -ea for java VM to enable assertions.
  */

 import java.lang.Long;
 import java.lang.Integer;
 import java.io.BufferedReader;
 import java.io.InputStreamReader;
 import java.io.FileInputStream;
 import java.io.BufferedWriter;
 import java.io.OutputStreamWriter;
 import java.io.FileOutputStream;
 import java.io.PrintWriter;
 import java.io.FileNotFoundException;
 import java.io.UnsupportedEncodingException;
 import java.util.StringTokenizer;
 import java.util.concurrent.*;

 abstract class SparseMatrixQ4 {
     int num_vertices; // Number of vertices in the graph
     int num_edges;    // Number of edges in the graph
     int [] source;    // sources from graphs
     int [] destination; // destinations from graph

     // Auxiliary in preparation of PageRank iteration: pre-calculate the
     // out-degree (number of outgoing edges) for each vertex
     abstract void calculateOutDegree(int outdeg[]);

     // Perform one PageRank iteration.
     //    a: damping factor
     //    in[]: previous PageRank values, read-only
     //    out[]: new PageRank values, initialised to zero
     //    outdeg[]: values pre-calculated by calculateOutDegree()
     //    start: value for start of array
     //    end: value for end of array
     abstract void iterate(double a, double[] in, double[] out, int outdeg[], int start, int end);
 }

 // This class represents the adjacency matrix of a graph as a sparse matrix
// in coordinate format (COO)
 class SparseMatrixCOOQ4 extends SparseMatrixQ4 {
     SparseMatrixCOOQ4(String file) {
         try {
             InputStreamReader is
                     = new InputStreamReader(new FileInputStream(file), "UTF-8");
             BufferedReader rd = new BufferedReader(is);
             readFile(rd);
         } catch (FileNotFoundException e) {
             System.err.println("File not found: " + e);
             return;
         } catch (UnsupportedEncodingException e) {
             System.err.println("Unsupported encoding exception: " + e);
             return;
         } catch (Exception e) {
             System.err.println("Exception: " + e);
             return;
         }
     }

     int getNext(BufferedReader rd) throws Exception {
         String line = rd.readLine();
         if (line == null)
             throw new Exception("premature end of file");
         return Integer.parseInt(line);
     }

     void getNextPair(BufferedReader rd, int pair[]) throws Exception {
         String line = rd.readLine();
         if (line == null)
             throw new Exception("premature end of file");
         StringTokenizer st = new StringTokenizer(line);
         pair[0] = Integer.parseInt(st.nextToken());
         pair[1] = Integer.parseInt(st.nextToken());
     }

     void readFile(BufferedReader rd) throws Exception {
         String line = rd.readLine();
         if (line == null)
             throw new Exception("premature end of file");
         if (!line.equalsIgnoreCase("COO"))
             throw new Exception("file format error -- header");

         num_vertices = getNext(rd);
         num_edges = getNext(rd);

         source = new int[num_edges];
         destination = new int[num_edges];

         int edge[] = new int[2];
         for (int i = 0; i < num_edges; ++i) {
             getNextPair(rd, edge);
             source[i] = edge[0];
             destination[i] = edge[1];
         }
     }

     // Auxiliary function for PageRank calculation
     void calculateOutDegree(int outdeg[]) {
         for (int i = 0; i < num_edges; i++) {
             outdeg[source[i]]++;
         }
     }

     void iterate(double a, double[] in, double[] out, int outdeg[], int start, int end) {
         for (int i = 0; i < num_edges; i++) {
             if (outdeg[destination[i]] != 0) {
                 out[source[i]] += a * (in[destination[i]] / outdeg[destination[i]]);
             } else {
                 out[source[i]] += a * in[destination[i]];
             }
         }
     }
 }

 // This class represents the adjacency matrix of a graph as a sparse matrix
// in compressed sparse rows format (CSR), where a row index corresponds to
 class SparseMatrixCSRQ4 extends SparseMatrixQ4 {
     SparseMatrixCSRQ4(String file) {
         try {
             InputStreamReader is
                     = new InputStreamReader(new FileInputStream(file), "UTF-8");
             BufferedReader rd = new BufferedReader(is);
             readFile(rd);
         } catch (FileNotFoundException e) {
             System.err.println("File not found: " + e);
             return;
         } catch (UnsupportedEncodingException e) {
             System.err.println("Unsupported encoding exception: " + e);
             return;
         } catch (Exception e) {
             System.err.println("Exception: " + e);
             return;
         }
     }

     int getNext(BufferedReader rd) throws Exception {
         String line = rd.readLine();
         if (line == null)
             throw new Exception("premature end of file");
         return Integer.parseInt(line);
     }

     void readFile(BufferedReader rd) throws Exception {
         String line = rd.readLine();
         if (line == null)
             throw new Exception("premature end of file");
         if (!line.equalsIgnoreCase("CSR") && !line.equalsIgnoreCase("CSC-CSR"))
             throw new Exception("file format error -- header");

         num_vertices = getNext(rd);
         num_edges = getNext(rd);

         source = new int[num_vertices + 1];
         destination = new int[num_edges];
         int destination_count = 0;

         for (int i = 0; i < num_vertices; i++) {
             line = rd.readLine();
             if (line == null)
                 throw new Exception("premature end of file");
             String elm[] = line.split(" ");
             assert Integer.parseInt(elm[0]) == i : "Error in CSR file";
             if (elm.length == 1) {
                 source[i] = destination_count;
             }
             for (int j = 1; j < elm.length; ++j) {
                 int dst = Integer.parseInt(elm[j]);
                 // TODO:
                 //    Record an edge from source i to destination dst
                 destination[destination_count] = dst;
                 if (j == 1) {
                     source[i] = destination_count;
                 }
                 destination_count++;
             }
         }
         source[source.length - 1] = destination_count;
     }

     // Auxiliary function for PageRank calculation
     void calculateOutDegree(int outdeg[]) {
     }

     void iterate(double a, double[] in, double[] out, int outdeg_unused[], int start, int end) {
         int outdeg = 0;
         for (int i = 0; i < num_vertices; i++) {
             for (int j = source[i]; j < source[i + 1]; j++) {
                 outdeg = source[i + 1] - source[i];
                 out[destination[j]] += a * (in[i] / outdeg);

             }
         }
     }
 }

 // This class represents the adjacency matrix of a graph as a sparse matrix
// in compressed sparse columns format (CSC). The incoming edges for each
// vertex are listed.
 class SparseMatrixCSCQ4 extends SparseMatrixQ4 {
     SparseMatrixCSCQ4(String file) {
         try {
             InputStreamReader is
                     = new InputStreamReader(new FileInputStream(file), "UTF-8");
             BufferedReader rd = new BufferedReader(is);
             readFile(rd);
         } catch (FileNotFoundException e) {
             System.err.println("File not found: " + e);
             return;
         } catch (UnsupportedEncodingException e) {
             System.err.println("Unsupported encoding exception: " + e);
             return;
         } catch (Exception e) {
             System.err.println("Exception: " + e);
             return;
         }
     }

     public int[] getSource() {
         return source;
     }

     int getNext(BufferedReader rd) throws Exception {
         String line = rd.readLine();
         if (line == null)
             throw new Exception("premature end of file");
         return Integer.parseInt(line);
     }

     void readFile(BufferedReader rd) throws Exception {
         String line = rd.readLine();
         if (line == null)
             throw new Exception("premature end of file");
         if (!line.equalsIgnoreCase("CSC") && !line.equalsIgnoreCase("CSC-CSR"))
             throw new Exception("file format error -- header");

         num_vertices = getNext(rd);
         num_edges = getNext(rd);

         // TODO: allocate data structures
         destination = new int[num_vertices + 1];
         source = new int[num_edges];
         int source_count = 0;

         for (int i = 0; i < num_vertices; ++i) {
             line = rd.readLine();
             if (line == null)
                 throw new Exception("premature end of file");
             String elm[] = line.split(" ");
             assert Integer.parseInt(elm[0]) == i : "Error in CSC file";
             if (elm.length == 1) {
                 destination[i] = source_count;
             }
             for (int j = 1; j < elm.length; ++j) {
                 int src = Integer.parseInt(elm[j]);
                 // TODO:
                 //    Record an edge from source src to destination
                 source[source_count] = src;
                 if (j == 1) {
                     destination[i] = source_count;
                 }
                 source_count++;
             }
         }
         destination[destination.length - 1] = source_count;
     }

     // Auxiliary function for PageRank calculation
     void calculateOutDegree(int outdeg[]) {
         // TODO:
         //    Calculate the out-degree for every vertex, i.e., the
         //    number of edges where a vertex appears as a source vertex.
         for (int i = 0; i < num_edges; i++) {
             ++outdeg[source[i]];
         }

     }

     void iterate(double a, double[] in, double[] out, int outdeg[], int start, int end) {
         // TODO:
         //    Iterate over all edges in the sparse matrix and calculate
         //    the contribution to the new PageRank value of a destination
         //    vertex made by the corresponding source vertex
         double [] contribution = new double[num_vertices];
         for(int i = 0; i < num_vertices; i++)
             contribution[i] = a * (in[i]/outdeg[i]);

         for (int i = start; i < end; i++) {
             for (int j = destination[i]; j < destination[i + 1]; j++) {
                 out[i] += contribution[source[j]];
             }
         }
     }
 }

 class PageRankQ4 implements Runnable {
     // Final variables
     static final double d = 0.85;
     static final double tol = 1e-7;
     static final int max_iter = 100;
     static final boolean verbose = true;

     // Variables for page rank
     static int n;
     static double x[];
     static double v[];
     static double y[];
     static int outdeg[];
     static double delta = 2;
     static int iter = 0;

     // Variables for time
     static long tm_start;

     // Variables relating to matrix
     static String outputFile;
     static SparseMatrixQ4 matrix;
     static int numVerts;
     static int numEdges;

     // Concurrent variables
     static boolean flag = false;
     static CyclicBarrier barrier;
     static int num_threads;
     static int [] start;
     static int [] end;

     @Override
     public void run() {
         while (iter < max_iter && delta > tol) {
             try {
                 barrier.await();
             } catch (InterruptedException ex) {
                 return;
             } catch (BrokenBarrierException ex) {
                 return;
             }
             if (flag) {
                 break;
             }
             // Power iteration step.
             // 1. Transferring weight over out-going links (summation part)
             int curThread = Integer.parseInt(Thread.currentThread().getName());
             //matrix.iterate(d, x, y, outdeg, curThread * ((n-1) / num_threads), (curThread + 1) * ((n-1) / num_threads));
             matrix.iterate(d, x, y, outdeg, start[curThread], end[curThread]);

             // 2. Constants (1-d)v[i] added in separately.
             try {
                 barrier.await();
             } catch (InterruptedException ex) {
                 return;
             } catch (BrokenBarrierException ex) {
                 return;
             }

             if (Thread.currentThread().getName().equals("0")) {
                 double w = 1.0 - sum(y, n); // ensure y[] will sum to 1
                 for (int i = 0; i < n; ++i)
                     y[i] += w * v[i];

                 // Calculate residual error
                 delta = normdiff(x, y, n);
                 iter++;

                 // Rescale to unit length and swap x[] and y[]
                 w = 1.0 / sum(y, n);
                 for (int i = 0; i < n; ++i) {
                     x[i] = y[i] * w;
                     y[i] = 0.;
                 }

                 if (delta < tol)
                     flag = true;

                 double tm_step = (double) (System.nanoTime() - tm_start) * 1e-9;
                 if (verbose)
                     System.err.println("iteration " + iter + ": delta=" + delta
                             + " xnorm=" + sum(x, n)
                             + " time=" + tm_step + " seconds");
                 tm_start = System.nanoTime();
             }
         }

         if (delta > tol)
             System.err.println("Error: solution has not converged.");

         if (Thread.currentThread().getName().equals("0")) {
             try {
                 barrier.await();
             } catch (InterruptedException ex) {
                 return;
             } catch (BrokenBarrierException ex) {
                 return;
             }
         }
         // Dump PageRank values to file
         writeToFile(outputFile, x, n);
     }

     public static void main(String args[]) {
         if (args.length < 5) {
             System.err.println("Usage: java pagerank graph.coo graph.csr graph.csc threads outputfile");
             return;
         }

         tm_start = System.nanoTime();

         System.out.println("COO: " + args[0]);
         System.out.println("CSR: " + args[1]);
         System.out.println("CSC: " + args[2]);

         // SparseMatrixQ3 matrix = new SparseMatrixCOOQ3( args[0] );
         // SparseMatrixQ3 matrix = new SparseMatrixCSRQ3( args[1] );
         matrix = new SparseMatrixCSCQ4(args[2]);

         num_threads = Integer.parseInt(args[3]);
         System.out.println("Number of threads: " + num_threads);
         outputFile = args[4];

         long tm_end = System.nanoTime();
         double tm_input = (double) (tm_end - tm_start) * 1e-9;
         tm_start = tm_end;
         System.out.println("Reading input: " + args[0] + " seconds");

         n = matrix.num_vertices;
         x = new double[n];
         v = new double[n];
         y = new double[n];
         barrier = new CyclicBarrier(num_threads);

         for (int i = 0; i < n; ++i) {
             x[i] = v[i] = ((double) 1) / (double) n;
             y[i] = 0;
         }

         outdeg = new int[n];
         matrix.calculateOutDegree(outdeg);

         //Pre calc number of edges to calc,

         start = new int [num_threads];
         end = new int[num_threads];

         int totEdgesToCalcPerThread = matrix.num_edges/num_threads;
         int startCount = 0;
         int endCount = 0;

         for(int i=0; i<num_threads; i++){
             if(i == num_threads-1){
                 start[i] = startCount;
                 end[i] = matrix.destination.length -1;
                 break;
             }
             for(int j=startCount; j<matrix.destination.length; j++){
                 if(matrix.destination[endCount] - matrix.destination[startCount] >= totEdgesToCalcPerThread)
                 {
                     start[i] = startCount;
                     end[i] = endCount;
                     startCount = endCount+1;
                     break;
                 }
                 else{
                     endCount++;
                 }
             }
         }

         for(int i=0; i<start.length; i++)
             System.err.println("Start source:\t" + start[i] + "\nEnd source:\t" + end[i]);


         double tm_init = (double) (System.nanoTime() - tm_start) * 1e-9;
         System.err.println("Initialisation: " + tm_init + " seconds");
         tm_start = System.nanoTime();

         // Start Threads
         Thread threads[] = new Thread[num_threads];
         for (int i = 0; i < num_threads; i++) {
             threads[i] = new Thread(new PageRankQ4());
             threads[i].setName(Integer.toString(i));
             threads[i].start();
         }
     }

     static double sum(double[] a, int n) {
         double d = 0.;
         double err = 0.;
         for (int i = 0; i < n; ++i) {
             // The code below achieves
             // d += a[i];
             // but does so with high accuracy
             double tmp = d;
             double y = a[i] + err;
             d = tmp + y;
             err = tmp - d;
             err += y;

         }
         return d;
     }

     static double normdiff(double[] a, double[] b, int n) {
         double d = 0.;
         double err = 0.;
         for (int i = 0; i < n; ++i) {
             // The code below achieves
             // d += Math.abs(b[i] - a[i]);
             // but does so with high accuracy
             double tmp = d;
             double y = Math.abs(b[i] - a[i]) + err;
             d = tmp + y;
             err = tmp - d;
             err += y;
         }
         return d;
     }

     static void writeToFile(String file, double[] v, int n) {
         try {
             OutputStreamWriter os
                     = new OutputStreamWriter(new FileOutputStream(file), "UTF-8");
             BufferedWriter wr = new BufferedWriter(os);
             writeToBuffer(wr, v, n);
         } catch (FileNotFoundException e) {
             System.err.println("File not found: " + e);
             return;
         } catch (UnsupportedEncodingException e) {
             System.err.println("Unsupported encoding exception: " + e);
             return;
         } catch (Exception e) {
             System.err.println("Exception: " + e);
             return;
         }
     }

     static void writeToBuffer(BufferedWriter buf, double[] v, int n) {
         PrintWriter out = new PrintWriter(buf);
         for (int i = 0; i < n; ++i)
             out.println(i + " " + v[i]);
         out.close();

     }
 }
