package Location::GeoTool;

################################################################
#
#  Geometric Functions 
#  Location::GeoTool
#  

use 5.008;
use strict;
use warnings;
use vars qw($VERSION @ISA @EXPORT);
$VERSION = 1.01;

use Math::Trig qw(asin tan);

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(
	coordformat
	dir_dist2point
	point2dir_dist
	datumchange
	LG_WGS84
	LG_TOKYO
	LG_FMT
	LG_FMT_MAPION
	LG_FMT_DMSN
	LG_FMT_SEC
	LG_FMT_DEG
	LG_FMT_RAD
	LG_FMT_GPSONE
);

my $pi  = 4 * atan2(1,1); 								# PI
my $rd  = $pi / 180;      								# [radian/degree]

################################################################
# Constant Methods             #
################################

sub LG_WGS84{0};
sub LG_TOKYO{1};
sub LG_FMT{0};
sub LG_FMT_MAPION{1};
sub LG_FMT_DMSN{2};
sub LG_FMT_SEC{3};
sub LG_FMT_DEG{4};
sub LG_FMT_RAD{5};
sub LG_FMT_GPSONE{6};

################################################################
# External Methods             #
################################

################################
# Change the format of coordinate

sub coordformat
{
	my ($coord,$farg,$targ,$pon) = @_[0..3];
	my $s = ($coord =~ /^(\+|-)/) ? $1 : "+";
	my ($dd, $mm, $ss, $nn);
	if ($farg == LG_FMT)
	{
		$coord =~ /(\d{3})(\d{2})(\d{2})(\d{3})$/;
		$dd = $1;
		$mm = $2;
		$ss = $3;
		$nn = $4;
	} 
	elsif ($farg == LG_FMT_MAPION) 
	{
		$coord =~ /(\d+)\/(\d+)\/(\d+)\.(\d+)$/;
		$dd = $1;
		$mm = $2;
		$ss = $3;
		$nn = $4;
	} 
	elsif ($farg == LG_FMT_DMSN)
	{
		$coord =~ /(\d+)(\d{2})(\d{2})\.(\d+)$/;
		$dd = $1;
		$mm = $2;
		$ss = $3;
		$nn = $4;
	} 
	elsif (($farg == LG_FMT_SEC) || ($farg == LG_FMT_DEG) || ($farg == LG_FMT_RAD))
	{
		if ($farg == LG_FMT_SEC) 
		{
			$nn = $coord/3600;
		} 
		elsif ($farg == LG_FMT_DEG) 
		{
			$nn = $coord;
		} 
		else 
		{
			$nn = $coord / $rd;
		}
		$dd = int($nn);
		$nn = ($nn-$dd) * 60;
		$mm = int($nn);
		$nn = ($nn-$mm) * 60;
		$ss = int($nn);
		$nn = int(($nn-$ss) * 1000);
		$nn = sprintf("%03d",$nn);
	} 
	elsif ($farg == LG_FMT_GPSONE) 
	{
		$coord =~ /(\d+)\.(\d+)\.(\d+)\.(\d+)$/;
		$dd = $1;
		$mm = $2;
		$ss = $3;
		$nn = $4;
	}

	$mm = '00'.$mm;
	$ss = '00'.$ss;
	$nn = $nn.'000';
	$mm = substr($mm,length($mm)-2);
	$ss = substr($ss,length($ss)-2);
	$nn = substr($nn,0,3);

	my $ret;
	if ($targ == LG_FMT)
	{
		$pon = 1;
		$ret = sprintf("%03d%02d%02d%03d",$dd,$mm,$ss,$nn);
	} 
	elsif ($targ == LG_FMT_MAPION) 
	{
		$ret = "$dd/$mm/$ss.$nn";
	} 
	elsif ($targ == LG_FMT_DMSN)
	{
		$ret = "$dd$mm$ss.$nn";
	} 
	elsif (($targ == LG_FMT_SEC) || ($targ == LG_FMT_DEG) || ($targ == LG_FMT_RAD))
	{
		$ss += ($dd * 3600 + $mm * 60);
		$ret = "$ss.$nn";
		if ($targ == LG_FMT_DEG) 
		{
			$ret = $ret/3600;
		} 
		if ($targ == LG_FMT_RAD) 
		{
			$ret = $ret*$rd/3600;
		} 
	} 
	elsif ($targ == LG_FMT_GPSONE)
	{
		$ret = "$dd.$mm.$ss.$nn";
	}
	
	# Set Plus/Minus sign if wanted
	if (($pon) && ($pon == 1)) 
	{
		$ret = $s.$ret;
	}

	return $ret;
}

################################
# Calcurate a point that has the direction 
# and the distance from other point

sub dir_dist2point
{
	my ($lat,$lon) = map {coordformat($_,LG_FMT,LG_FMT_RAD)} @_[0..1];
	my ($dir,$dis,$datum) = @_[2..4];					# Direction,Distance,Datum

	$dir = $dir * $rd;												# Change to radian

	my $ellip = ellipsoid($datum);
	my $a = $ellip->{'a'};										# Equatorial Radius
	my $f = $ellip->{'f'};										# Flattening
	my $r = 1 - $f;
	my $tu = $r * tan($lat);
	my $sf = sin($dir);
	my $cf = cos($dir);
	my $b = ($cf == 0) ? 0.0 : 2.0 * atan2($tu,$cf);

	my $cu = 1.0 / sqrt(1 + $tu**2);
	my $su = $tu * $cu;
	my $sa = $cu * $sf;
	my $c2a = 1 - $sa**2;
	my $x = 1.0 + sqrt(1.0 + $c2a * (1.0/($r**2)-1.0));
	$x = ($x - 2.0) / $x;

	my $c = 1.0 - $x;
	$c = ($x**2 / 4.0 + 1.0) / $c;
	my $d = (0.375 * $x**2 - 1.0)* $x;
	$tu = $dis / ($r * $a * $c);
	my $y = $tu;
	$c = $y + 1;

	my ($sy,$cy,$cz,$e) = ();
	while (abs($y - $c) > 0.00000000005)
	{
		$sy = sin($y);
		$cy = cos($y);
		$cz = cos($b + $y);
		$e = 2.0 * $cz**2 -1.0;
		$c = $y;
		$x = $e * $cy;
		$y = $e + $e - 1;
		$y = ((($sy**2 * 4.0 - 3.0) * $y * $cz * $d / 6.0 + $x) * $d / 4.0 - $cz) * $sy * $d + $tu;
	}
		
	$b = $cu * $cy * $cf - $su * $sy;
	$c = $r * sqrt($sa**2 + $b**2);
	$d = $su * $cy + $cu * $sy * $cf;
	my $rlat = atan2($d,$c);

	$c = $cu * $cy - $su * $sy * $cf;
	$x = atan2($sy * $sf, $c); 
	$c = ((-3.0 * $c2a + 4.0) * $f + 4.0) * $c2a * $f / 16.0;
	$d = (($e * $cy * $c + $cz) * $sy * $c + $y) * $sa;
	my $rlon = $lon + $x - (1.0 - $c) * $d * $f;

	return map {coordformat($_,LG_FMT_RAD,LG_FMT,1)} ($rlat,$rlon);
}

################################
# Calcurate distance and direction of points

sub point2dir_dist
{
	my ($lat,$lon,$tlat,$tlon) = map {coordformat($_,LG_FMT,LG_FMT_RAD)} @_[0..3];
	my $datum = $_[4];

	return (180,0) if (($lat == $tlat) && ($lon == $tlon));

	my $ellip = ellipsoid($datum);
	my $a = $ellip->{'a'};										# Equatorial Radius
	my $f = $ellip->{'f'};										# Flattening
  my $e2  = 2*$f - $f*$f;   								# Square of Eccentricity

	my $r = 1 - $f;
	my $tu1 = $r * tan($lat);
	my $tu2 = $r * tan($tlat);
	my $cu1 = 1.0 / sqrt(1.0 + $tu1**2);
	my $su1 = $cu1 * $tu1;
	my $cu2 = 1.0 / sqrt(1.0 + $tu2**2); 
	my $s1 = $cu1 * $cu2;
	my $b1 = $s1 * $tu2;
	my $f1 = $b1 * $tu1;
	my $x = $tlon - $lon;
	my $d = $x + 1;														# Force one pass
	my $iter =1;
	my ($sx,$cx,$sy,$cy,$y,$sa,$c2a,$cz,$e,$c)=();
	while ((abs($d - $x) > 0.00000000005) && ($iter < 100))
	{
		$iter++;
		$sx = sin($x);
		$cx = cos($x);
		$tu1 = $cu2 * $sx;
		$tu2 = $b1 - $su1 * $cu2 * $cx;
		$sy = sqrt($tu1**2 + $tu2**2);
		$cy = $s1 * $cx + $f1;
		$y = atan2($sy,$cy);
		$sa = $s1 * $sx / $sy;
		$c2a = 1 - $sa**2;
		$cz = $f1 + $f1;
		if ($c2a > 0.0)
		{
			$cz = $cy - $cz / $c2a;
		}
		$e = $cz**2 * 2.0 - 1.0;
		$c = ((-3.0 * $c2a + 4.0) * $f + 4.0) * $c2a * $f / 16.0;
		$d = $x;
		$x = (($e * $cy * $c + $cz) * $sy * $c + $y) * $sa;
		$x = (1.0 - $c) * $x * $f + $tlon - $lon;
	}
	my $dir = atan2($tu1,$tu2) / $rd;
	$x = sqrt((1 / ($r**2) -1) * $c2a +1);
	$x += 1;
	$x = ($x - 2.0) / $x;
	$c = 1.0 - $x;
	$c = ($x**2 / 4.0 + 1.0) / $c;
	$d = (0.375 * $x**2 - 1.0) * $x;
	$x = $e * $cy;
	my $dis = (((($sy**2 * 4.0 - 3.0) * (1.0 - $e - $e) * $cz * $d / 6.0 - $x) * $d / 4.0 + $cz) * $sy * $d + $y) * $c * $a * $r;
	
	return ($dir,$dis);
}

################################
# Change the coordinate to different datum

sub datumchange{
	my ($lat,$lon,$h,$from,$to) = @_;
	my ($b,$l) = map { coordformat($_,LG_FMT,LG_FMT_DEG) } ($lat,$lon);
	$h = eval($h)/100;

	my $fellip = ellipsoid($from);
	my $tellip = ellipsoid($to);
	my @a = ($fellip->{'a'},$tellip->{'a'});			# Equatorial Radius
	my @f = ($fellip->{'f'},$tellip->{'f'});			# Flattening

	my $dx = $tellip->{'dx'}-$fellip->{'dx'};
	my $dy = $tellip->{'dy'}-$fellip->{'dy'};
	my $dz = $tellip->{'dz'}-$fellip->{'dz'};

	($b, $l, $h) = molodensky($b, $l, $h, $a[0], $f[0], $a[1], $f[1],$dx,$dy,$dz);
	my ($tlat,$tlon) = map {coordformat($_,LG_FMT_DEG,LG_FMT,1)} ($b,$l);
	$h = sprintf("%+06d",int($h*100));
	return ($tlat,$tlon,$h);
}

################################################################
# Internal Methods             #
################################

################################
# Constants of each datum

sub ellipsoid
{
	my $datum = $_[0] ? $_[0] : LG_WGS84;
	my @ellip;
	$ellip[LG_WGS84] = [6378137,(1 / 298.257223),0,0,0];
	$ellip[LG_TOKYO] = [6377397.155,(1 / 299.152813),148,-507,-681];

	return 
	{
		'a' => $ellip[$datum]->[0],						# Equatorial Radius
		'f' => $ellip[$datum]->[1],						# Flattening
		'dx' => $ellip[$datum]->[2],					# delta x parameter
		'dy' => $ellip[$datum]->[3],					# delta y parameter
		'dz' => $ellip[$datum]->[4]						# delta z parameter
	};
}

################################
# Main procedure of Molodensky method

sub molodensky {
	my($b, $l, $h, $a, $f, $a_, $f_,$dx,$dy,$dz) = @_;
	my($bda, $e2, $da, $df, $db, $dl, $dh);
	my($sb, $cb, $sl, $cl, $rn, $rm);

	$b *= $rd;
	$l *= $rd;

	$e2 = 2*$f - $f*$f; 						# Square of Eccentricity
	$bda = 1- $f;       						# Polar Radius / Equatorial Radius
	($da, $df) = ($a_-$a, $f_-$f);
	($sb, $cb, $sl, $cl) = (sin($b), cos($b), sin($l), cos($l));

	$rn = 1 / sqrt(1 - $e2*$sb*$sb); 
	$rm = $a * (1 - $e2) * $rn * $rn * $rn;
	$rn *= $a;

	# Calcurating Delta Value
  $db = -$dx*$sb*$cl - $dy*$sb*$sl + $dz*$cb
		+ $da*$rn*$e2*$sb*$cb/$a + $df*($rm/$bda+$rn*$bda)*$sb*$cb;
	$db /= $rm + $h;
	$dl = -$dx*$sl + $dy*$cl;
	$dl /= ($rn+$h) * $cb;
	$dh = $dx*$cb*$cl + $dy*$cb*$sl + $dz*$sb
		- $da*$a/$rn + $df*$bda*$rn*$sb*$sb;

	return (($b+$db)/$rd, ($l+$dl)/$rd, $h+$dh);
}

1;
__END__

=head1 NAME

Location::GeoTool - Perl extension for Geometry processing

=head1 SYNOPSIS

  use Location::GeoTool;
  
  my $ret = coordformat($coord,$farg,$targ,$pon);
  my ($tlat,$tlon) = dir_dist2point($lat,$lon,$dir,$dis,$datum);
  my ($dir,$dis) = point2dir_dist($lat,$lon,$tlat,$tlon,$datum);
  my ($tlat,$tlon) = datumchange($lat,$lon,$h,$from,$to);

=head1 DESCRIPTION

=head2 EXPORT

  coordformat
  dir_dist2point
  point2dir_dist
  datumchange

=head3 coordformat

  Change the format of coordinate.
  Standard format is pdddmmssnnn.
  (p: plus/minus, d:degree, m:minute, s:second n: micro second)

  Usage:
    $ret = coordformat($coord,$farg,$targ,$pon);
  Arguments:
    $coord : Original Formatted Coordinate
    $farg  : Original Format ID
    $targ  : Wanted Format ID
    $pon   : On/Off of plus/minus
  Return Values:
    $ret   : Wanted Formatted Coordinate
  Format ID:
    LG_FMT         : Standard Format (pdddmmssnnn)
    LG_FMT_MAPION  : Mapion URL Format (ddd/mm/ss.nnn)
    LG_FMT_DMSN    : dddmmss.nnn
    LG_FMT_SEC     : Second (ssssss.sss)
    LG_FMT_DEG     : Degree (ddd.dddddddd)
    LG_FMT_RAD     : Radian
    LG_FMT_GPSONE  : ddd.mm.ss.nnn

=head3 dir_dist2point 

  Calcurate a point that has the direction and the distance from 
  other point.
  Original Version is in JavaScript,
  http://williams.best.vwh.net/gccalc.htm
  coding by Mr. Ed Williams.

  Coordinates format is : pdddmmssnnn

  Usage:
    ($tlat,$tlon) = dir_dist2point($lat,$lon,$dir,$dis,$datum);
  Arguments:
    $lat   : Latitude
    $lon   : Longitude
    $dir   : Direction (degree, North and East are Positive.)
    $dis   : Distance (m)
    $datum : LG_WGS84, LG_TOKYO
  Return Values:
    $tlat  : Latitude of calcurated position
    $tlon  : Longitude of calcurated position

=head3 point2dir_dist

  Give two points by coordinates, then calcurate distance 
  and direction of them.
  Original Version is in JavaScript,
  http://williams.best.vwh.net/gccalc.htm
  coding by Mr. Ed Williams.

  Coordinates format is : pdddmmssnnn

  Usage:
    ($dir,$dis) = point2dir_dist($lat,$lon,$tlat,$tlon,$datum);
  Arguments:
    $lat   : Latitude of From Point
    $lon   : Longitude of From Point
    $tlat  : Latitude of To Point
    $tlon  : Longitude of To Point
    $datum : LG_WGS84, LG_TOKYO
  Return Values:
    $dir   : Direction (degree, North and East are Positive.)
    $dis   : Distance (m)

=head3 datumchange

  Change the coordinate to different datum.
  Logic is based on Molodensky method. 
  Original Version is Perl CGI
  http://member.nifty.ne.jp/Nowral/02_DATUM/Molodensky.html
  coding by Nowral.

  Coordinates format is : pdddmmssnnn

  Usage:
    ($tlat,$tlon,$th) = datumchange($lat,$lon,$h,$from,$to);
  Arguments:
    $lat  : Latitude
    $lon  : Longitude
    $h¡¡  : Altitude (cm)
    $from : Original Datum
    $to   : Wanted Datum
  Return Values:
    $tlat : Latitude in wanted Datum
    $tlon : Longitude in wanted Datum
    $th¡¡ : Altitude in wanted Datum

=head1 DEPENDENCIES

Math::Trig

=head1 SEE ALSO

dir_dist2point, dir_dist2point function is based on javascript program
could be seen in
http://williams.best.vwh.net/gccalc.htm

datumchange and molodensly function is based on perl program could 
be seen in
http://member.nifty.ne.jp/Nowral/02_DATUM/Molodensky.html

Thanks for these site.

Support this module in SpaceTag Inc. web site : http://www.spacetag.jp/

=head1 AUTHOR

OHTSUKA Ko-hei, E<lt>kotsuka@spacetag.jpE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2003 by SpaceTag INC.,

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.1 or,
at your option, any later version of Perl 5 you may have available.


=cut
