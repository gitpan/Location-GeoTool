use strict;
use Test::More tests => 33;

BEGIN { use_ok 'Location::GeoTool' }

my %testcase = (
  61 => 
  {
    4 => ['��','E'],
    8 => ['����','NE'],
    16 => ['������','ENE'],
    32 => ['��������','NEbE']
  },
  124 => 
  {
    4 => ['��','E'],
    8 => ['����','SE'],
    16 => ['����','SE'],
    32 => ['��������','SEbE']
  },
  247 => 
  {
    4 => ['��','W'],
    8 => ['����','SW'],
    16 => ['������','WSW'],
    32 => ['������','WSW']
  },
  324 => 
  {
    4 => ['��','N'],
    8 => ['����','NW'],
    16 => ['����','NW'],
    32 => ['��������','NWbW']
  }
);

my $obj = Location::GeoTool->create_coord(35.12345,139.12345,'wgs84','degree');

foreach my $degree (keys %testcase)
{
  my $dir = $obj->direction_vector($degree,100);
  foreach my $mother (keys %{$testcase{$degree}})
  {
    my @string = @{$testcase{$degree}->{$mother}};
    is $dir->dir_string($mother,'jp'),$string[0];
    is $dir->dir_string($mother,'en'),$string[1];
  }
}

