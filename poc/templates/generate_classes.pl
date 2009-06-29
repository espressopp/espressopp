#!/usr/bin/perl

for ($i = 0; $i < 10; $i++) {

print <<EOT;
class Class_$i {
public:
  double get() {
      double sum = 0.0;
EOT

for ($j=0; $j < 3; $j++) {
    $r = int(rand(9));
    if ($r == 0) {
	print "sum += sin(sum);\n";
    } elsif ($r == 1) {
	print "sum += cos(sum);\n";
    } elsif ($r == 2) {
	print "sum += exp(sum);\n";
    } elsif ($r == 3) {
	print "sum += sqrt(sum);\n";
    } elsif ($r == 4) {
	print "sum += tan(sum);\n";
    } elsif ($r == 5) {
	print "sum += cosh(sum);\n";
    } elsif ($r == 6) {
	print "sum += sinh(sum);\n";
    } elsif ($r == 7) {
	print "sum += sum < 0.5 ? 0.0 : 1.0;\n";
    } elsif ($r == 8) {
	print "	for (int i = 0; i < 3; i++) { sum += int(sum); }\n";
    }
}

 print <<EOT;
    return sum;
  }
};
EOT

}
