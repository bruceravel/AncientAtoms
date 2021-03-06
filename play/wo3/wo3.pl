#!/usr/bin/perl -w

@sites = (["W", "w(1)", 0.2465, 0.0269, 0.2859],
	  ["W", "w(2)", 0.2538, 0.0353, 0.7808],
	  ["O", "ox1",  0.0025, 0.0350, 0.2106],
	  ["O", "ox2",  0.9974, 0.4636, 0.2161],
	  ["O", "oy1",  0.2840, 0.2605, 0.2848],
	  ["O", "oy2",  0.2099, 0.2586, 0.7318],
	  ["O", "oz1",  0.2827, 0.0383, 0.0046],
	  ["O", "oz2",  0.2856, 0.4840, 0.9944]);


foreach $site (@sites) {

  ($x,$y,$z) = ($site->[2], $site->[3], $site->[4]);
  printf "  %s\t%7.4f\t%7.4f\t%7.4f\t\t%s_1\n",
  $site->[0],  $x,  $y,  $z,  $site->[1];

  ($x,$y,$z) = (-$site->[2]+0.5, $site->[3]+0.5, -$site->[4]+0.5);
  printf "  %s\t%7.4f\t%7.4f\t%7.4f\t\t%s_2\n",
  $site->[0],  $x,  $y,  $z,  $site->[1];

  ($x,$y,$z) = (-$site->[2], -$site->[3], -$site->[4]);
  printf "  %s\t%7.4f\t%7.4f\t%7.4f\t\t%s_3\n",
  $site->[0],  $x,  $y,  $z,  $site->[1];

  ($x,$y,$z) = ($site->[2]+0.5, -$site->[3]+0.5, $site->[4]+0.5);
  printf "  %s\t%7.4f\t%7.4f\t%7.4f\t\t%s_4\n",
  $site->[0],  $x,  $y,  $z,  $site->[1];

};
