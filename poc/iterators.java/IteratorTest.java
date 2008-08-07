public class IteratorTest {
    private final static int N = 10000;
    private final static int NUM_TESTS = 10000;

     public static void main(String[] args) {
         System.out.println("Testing iterators...");

	 double[][] a = new double[N][3];
	 java.util.Random rand = new java.util.Random();

	 for (int i = 0; i < N; i++)
		 for (int j = 0; j < 3; j++)
			 a[i][j] = rand.nextDouble();
	 	
	 double[] sum = { 0.0, 0.0, 0.0 };
	 java.util.Date start = new java.util.Date();
	 for (int i= 0; i < NUM_TESTS; i++) {
	     for (double[] v:a) {
	    	 sum[0] += v[0];
	    	 sum[1] += v[1];
	    	 sum[2] += v[2];
	     }
	 }
	 java.util.Date end = new java.util.Date();
	 
	 double time = (end.getTime() - start.getTime()) / 1000.0;
	 
	 System.out.println("sum={" + 
			 			 sum[0] + ", " + 
			 			 sum[1] + ", " + 
			 			 sum[2] + "}");
	 System.out.println("time=" + time);

     }
 }
