use strict;
use warnings;

sub mer2AA {
my $str = split(/ /,$str);
my $str = shift;
$str =~ tr/0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19/A C D E F G H I K L M N P Q R S T V W Y/;
return $str
}

sub rot13 {
my $str = shift;
$str =~ tr/A-Za-z/N-ZA-Mn-za-m/;
return $str;
}

sub rot47 {
my $str = shift;
$str =~ tr/!-~/P-~!-O/;
return $str;
}

my $string = '1 2 3 4 5 6 19';

print "$string\n";
print mer2AA($string)."\n";
#print rot13($string)."\n";
#print rot47($string)."\n";
