#!/usr/bin/env perl

$| = 1; # autoflush

$SLEEP = 3600;

$i = 1;
$last = time();
while(1) {
  print "STARTING pollsim ", $i++, " ... ";
  $start = time();
  system("./pollsim.m > log.txt") == 0 or die;
  #system("scp index.html yootles.com:y/pollsim");
  print "DONE in ", ss(time()-$start), " (and now sleeping ",ss($SLEEP),")\n";
  sleep($SLEEP);
}

# double-digit: takes number from 0-99, returns 2-char string eg "03" or "42".
sub dd { my($n) = @_;  return ($n<=9 && $n>=0 ? "0".$n : $n); }

# Seconds to str: takes number of seconds, returns a string like 1d02h03:04:05
sub ss { my($s) = @_;
  my($d,$h,$m);
  my $incl = "s";

  if ($s < 0) { return "-".ss(-$s); }

  $m = int($s/60);
  if ($m > 0) { $incl = "ms"; }
  $s %= 60;
  $h = int($m/60);
  if ($h > 0) { $incl = "hms"; }
  $m %= 60;
  $d = int($h/24);
  if ($d > 0) { $incl = "dhms"; }
  $h %= 24;

  return ($incl=~"d" ? "$d"."d" : "").
         ($incl=~"h" ? dd($h)."h" : "").
         ($incl=~"m" ? dd($m).":" : "").
         ($incl!~"m" ? $s : dd($s))."s";
}
