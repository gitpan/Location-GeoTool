use Test::More tests => 5;
use lib "./lib";
use Location::GeoTool;

# test 1: use test
use_ok('Location::GeoTool');

# test 2: coordformat test
my $coord1 = "+1350000000";
ok ((coordformat($coord1,LG_FMT,LG_FMT_MAPION) eq "135/00/00.000") &&
(coordformat($coord1,LG_FMT,LG_FMT_DMSN) eq "1350000.000") && 
(coordformat($coord1,LG_FMT,LG_FMT_SEC) eq "486000.000") &&
(coordformat($coord1,LG_FMT,LG_FMT_DEG) == 135) &&
(abs(coordformat($coord1,LG_FMT,LG_FMT_RAD) - 2.35619449) < 0.00000001) &&
(coordformat($coord1,LG_FMT,LG_FMT_GPSONE,1) eq "+135.00.00.000"));

# test 3: point2dir_dist test
my ($dir,$dist) = point2dir_dist("+0350000000","+1350000000","+0360000000","+1360000000");
ok ((abs($dir - 38.98534584) < 0.00000001) &&
(abs($dist - 143321.5781) < 0.0001));

# test 4: dir_dist2point test
my ($tlat,$tlon) = dir_dist2point("+0350000000","+1350000000",45,100);
($tlat,$tlon) = map { coordformat($_,LG_FMT,LG_FMT_DEG) } ($tlat,$tlon);
ok ((abs($tlat - 35.0006372222222) < 0.00000001) && 
(abs($tlon - 135.000774444444) < 0.0000001));

# test 5: datumchange test
my $th;
($tlat,$tlon,$th) = datumchange("+0350000000","+1350000000","+000000",LG_WGS84, LG_TOKYO);
($tlat,$tlon) = map { coordformat($_,LG_FMT,LG_FMT_DEG) } ($tlat,$tlon);
ok ((abs($tlat - 34.9968036111111) < 0.00000001) && 
(abs($tlon - 135.002780555556) < 0.0000001) && 
(abs($th + 5198) <10));

