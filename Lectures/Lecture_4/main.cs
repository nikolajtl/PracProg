using static System.Console;
class main{
    static int Main(){
        double x=0;
        double y=x;
        x=9;
        int n=3;
        double[] a = new double[n];
        for (int i=0; i<a.Length; i++) a[i]=i+1;
        double[] b = a;
        b[0]=100;
        WriteLine($"a[0]={a[0]}");
        return 0;
    }
}