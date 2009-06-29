#!/usr/bin/perl

for ($i = 0; $i < 10; $i++) {
    for ($j = 0; $j < 10; $j++) {
	for ($k = 0; $k < 10; $k++) {
	    print("INSTANTIATE_CLASS($i,$j,$k);\n");
	}
    }
}
