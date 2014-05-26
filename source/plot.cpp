#include <iostream>
using namespace std;

#include "igpm_t2_lib.hpp"


////////////////////
#include "config.hpp"
#include "mesh.hpp"
#include "Tsolver.hpp"

//------------------------------------------------------------------------------
int main(int argc, char **argv) {

  if (argc<2) {  // argc
    cerr << "Please give name of config file as argument!" << endl;
    exit(1);
  }


  // open config file

  std::string configfilename(argv[1]);
  if (!singleton_config_file::instance().read_from_file(configfilename)) {
    cerr << "Unable to read from config file " << argv[1] << "!" << endl;
    exit(1);
  }

  cerr << "maxchildren = childdim = " << int(example_id_type::maxChildren) << endl;

  cerr << "------------------------------------------------------------------------------" << endl;

  /*
   * read a configuration (text)file (./config/igpm_t2_grid_example.singleton_config_file::instance())
   * file content:
   *
   * [section]                        # comment
   * key1 = data1                     # single data
   * key2 = data2.1,data2.2,data2.3   # list of data, comma seperated
   * # etc.
   *
   * use this as a control file, for input data and tasks,
   * build a config.file library for test cases, experimental data, etc.
   *
   * main function:
   * bool getValue( section, key, variable-for-data, default-value )
   */

  Tsolver solver;

  {
    int  reserve_base, reserve_ante, reserve_leaf;

    singleton_config_file::instance().getValue("cells","reserve_base",reserve_base,1000);
    singleton_config_file::instance().getValue("cells","reserve_ante",reserve_ante,25000);
    singleton_config_file::instance().getValue("cells","reserve_leaf",reserve_leaf,25000);

    solver.grid.reserve(reserve_base, reserve_ante, reserve_leaf);
  }

  std::string gridfiletemplate, outputdir, output_prefix, msafile, lagfile;
  if(!singleton_config_file::instance().getValue("general","gridfiletemplate",gridfiletemplate,"")) {
    cerr << "Essential option [general] gridfiletemplate missing in config file!" << endl;
    exit(2);
  }
  if(!singleton_config_file::instance().getValue("general","outputdir",outputdir,"")) {
    cerr << "Essential option [general] outputdir missing in config file!" << endl;
    exit(2);
  }

  if(!singleton_config_file::instance().getValue("general","outputprefix",output_prefix,"")) {
    cerr << "Essential option [general] outputdir missing in config file!" << endl;
    exit(2);
  }

  if(!singleton_config_file::instance().getValue("general","msafile",msafile,"")) {
    cerr << "Essential option [general] msafile missing in config file!" << endl;
    exit(2);
  }
  if(!singleton_config_file::instance().getValue("general","lagfile",lagfile,"")) {
    cerr << "Essential option [general] lagfile missing in config file!" << endl;
    exit(2);
  }

  cerr << "gridfiletemplate = " << gridfiletemplate << endl;
  cerr << "outputdir = " << outputdir << endl;
  cerr << "outputprefix = " << output_prefix << endl;
  cerr << "msafile = " << msafile << endl;

  solver.read_easymesh_triangles (gridfiletemplate);
  solver.shape.read_sc_msa (msafile);

  cerr << "number of blocks: "     << solver.grid.countBlocks() << endl;
  cerr << "number of nodes: "      << solver.grid.countNodes() << endl;
  cerr << "number of base cells: " << solver.grid.baseCells().size() << endl;

  // now all base cells are given, have a look (tecplot -> load data file)
  {
    std::string fname(outputdir);
    fname+="/grid2D_basecells.dat";
    solver.grid.writeCells(fname.c_str(), solver.grid.baseCells());
  }

  // finalize the setting, i.e., the grid implements the cell connections
  solver.grid.finalize();

  //
  solver.shape.initialize_quadrature ();
  solver.shape.initialize_mass ();

  solver.initialize_plotter();

  solver.read_problem_parameters_GENERAL (outputdir, output_prefix);


#if (EQUATION == MONGE_AMPERE_EQ)
  cerr << "EQUATION == MONGE_AMPERE_EQ" << endl;
//  solver.time_stepping_MA();
  int stabsign;
  double gamma, refine_eps, coarsen_eps;

  int level;
  solver.read_problem_parameters_MA(stabsign, gamma, refine_eps, coarsen_eps, level);

  //check level
  if (solver.levelmax < level)
	  level = solver.levelmax;

  solver.update_baseGeometry();

  solver.set_leafcellmassmatrix();

  solver.refine(level);

  solver.initializeLeafCellData_MA();

  unsigned int dofs_DG = 384;

  // init dof offsets
  solver.assignViews_MA(dofs_DG);

  //solution calculated in matlab
/*
  Eigen::VectorXd solution(15);
  solution(0) = 0.6362698395;
  solution(1) = 0.7664851639;
  solution(2) = 1.2835952523;
  solution(3) = 1.6277595302;
  solution(4) = 1.8357295980;
  solution(5) = 2.4859622739;
  solution(6) = 0.5585150961;
  solution(7) = 0.9978704414;
  solution(8) = 0.5585150961;
  solution(9) = 1.9968084718;
  solution(10) = 1.5002159990;
  solution(11) = -0.0110552831;
  solution(12) = 0.1191597512;
  solution(13) = -0.0110553040;
  solution(14) = 0.1191597511;
*/

  Eigen::VectorXd solution(561);
  solution(0) = -0.1102696421;
  solution(1) = -0.1102696321;
  solution(2) = -0.0761768988;
  solution(3) = -0.1102696421;
  solution(4) = -0.0761768988;
  solution(5) = -0.0420841744;
  solution(6) = -0.1594435990;
  solution(7) = -0.1594436089;
  solution(8) = -0.1348566108;
  solution(9) = -0.1594436386;
  solution(10) = -0.1348566306;
  solution(11) = -0.1594436386;
  solution(12) = -0.1594436089;
  solution(13) = -0.1348566306;
  solution(14) = -0.1348566108;
  solution(15) = -0.1929056608;
  solution(16) = -0.1929056901;
  solution(17) = -0.1761746398;
  solution(18) = -0.1929057393;
  solution(19) = -0.1761746792;
  solution(20) = -0.2140540993;
  solution(21) = -0.2105803088;
  solution(22) = -0.2034799094;
  solution(23) = -0.2071065379;
  solution(24) = -0.2000061289;
  solution(25) = -0.2140540605;
  solution(26) = -0.2140540700;
  solution(27) = -0.2034798508;
  solution(28) = -0.2034798703;
  solution(29) = -0.1929057393;
  solution(30) = -0.1929056901;
  solution(31) = -0.1761746792;
  solution(32) = -0.1929056608;
  solution(33) = -0.1761746398;
  solution(34) = -0.2140540993;
  solution(35) = -0.2140540700;
  solution(36) = -0.2034798703;
  solution(37) = -0.2034798508;
  solution(38) = -0.2071065379;
  solution(39) = -0.2105803088;
  solution(40) = -0.2000061289;
  solution(41) = -0.2034799094;
  solution(42) = -0.1929056509;
  solution(43) = -0.1761746201;
  solution(44) = -0.1761746201;
  solution(45) = -0.2067585294;
  solution(46) = -0.1951066506;
  solution(47) = -0.2069325238;
  solution(48) = -0.1834547916;
  solution(49) = -0.1952806550;
  solution(50) = -0.1831068147;
  solution(51) = -0.1662425142;
  solution(52) = -0.1832807933;
  solution(53) = -0.1493782333;
  solution(54) = -0.1664165027;
  solution(55) = -0.2064105920;
  solution(56) = -0.1947586935;
  solution(57) = -0.2065845508;
  solution(58) = -0.1949326622;
  solution(59) = -0.1490302956;
  solution(60) = -0.1256439224;
  solution(61) = -0.1492042546;
  solution(62) = -0.1022575686;
  solution(63) = -0.1258178912;
  solution(64) = -0.1019096706;
  solution(65) = -0.0710161741;
  solution(66) = -0.1020836097;
  solution(67) = -0.0401226962;
  solution(68) = -0.0711901230;
  solution(69) = -0.1486824371;
  solution(70) = -0.1252960442;
  solution(71) = -0.1488563565;
  solution(72) = -0.1254699734;
  solution(73) = -0.2060627318;
  solution(74) = -0.1944108144;
  solution(75) = -0.2062366520;
  solution(76) = -0.1827589165;
  solution(77) = -0.1945847444;
  solution(78) = -0.1824110960;
  solution(79) = -0.1655467568;
  solution(80) = -0.1825849964;
  solution(81) = -0.1657206669;
  solution(82) = -0.1987672717;
  solution(83) = -0.1905891742;
  solution(84) = -0.2024149919;
  solution(85) = -0.1942369041;
  solution(86) = -0.1658945963;
  solution(87) = -0.1829328557;
  solution(88) = -0.1660685454;
  solution(89) = -0.1834547916;
  solution(90) = -0.1951066506;
  solution(91) = -0.1952806550;
  solution(92) = -0.2067585294;
  solution(93) = -0.2069325238;
  solution(94) = -0.1831068147;
  solution(95) = -0.1947586935;
  solution(96) = -0.1949326622;
  solution(97) = -0.2064105920;
  solution(98) = -0.2065845508;
  solution(99) = -0.1493782333;
  solution(100) = -0.1662425142;
  solution(101) = -0.1664165027;
  solution(102) = -0.1832807933;
  solution(103) = -0.1827589165;
  solution(104) = -0.1944108144;
  solution(105) = -0.1945847444;
  solution(106) = -0.2060627318;
  solution(107) = -0.2062366520;
  solution(108) = -0.1824110960;
  solution(109) = -0.1905891742;
  solution(110) = -0.1942369041;
  solution(111) = -0.2024149919;
  solution(112) = -0.1486824371;
  solution(113) = -0.1655467568;
  solution(114) = -0.1657206669;
  solution(115) = -0.1825849964;
  solution(116) = -0.1022575686;
  solution(117) = -0.1256439224;
  solution(118) = -0.1258178912;
  solution(119) = -0.1490302956;
  solution(120) = -0.1492042546;
  solution(121) = -0.1019096706;
  solution(122) = -0.1252960442;
  solution(123) = -0.1254699734;
  solution(124) = -0.1488563565;
  solution(125) = -0.0401226962;
  solution(126) = -0.0710161741;
  solution(127) = -0.0711901230;
  solution(128) = -0.1020836097;
  solution(129) = -0.1658945963;
  solution(130) = -0.1660685454;
  solution(131) = -0.1829328557;
  solution(132) = -0.2060627220;
  solution(133) = -0.2133582316;
  solution(134) = -0.2098844020;
  solution(135) = -0.2097104718;
  solution(136) = -0.2098844020;
  solution(137) = -0.2097104718;
  solution(138) = -0.2137061305;
  solution(139) = -0.2102323201;
  solution(140) = -0.2100583514;
  solution(141) = -0.2104063046;
  solution(142) = -0.2138800857;
  solution(143) = -0.2138801051;
  solution(144) = -0.2102323201;
  solution(145) = -0.2137061305;
  solution(146) = -0.2100583514;
  solution(147) = -0.2138801051;
  solution(148) = -0.2138800857;
  solution(149) = -0.2104063046;
  solution(150) = -0.2137061207;
  solution(151) = -0.2135321712;
  solution(152) = -0.2135321712;
  solution(153) = -0.0404706338;
  solution(154) = -0.0713640918;
  solution(155) = -0.0402966551;
  solution(156) = -0.1026055459;
  solution(157) = -0.1259918798;
  solution(158) = -0.1024315473;
  solution(159) = -0.0408186508;
  solution(160) = -0.0717120890;
  solution(161) = -0.0406446324;
  solution(162) = -0.0715380805;
  solution(163) = -0.1497262499;
  solution(164) = -0.1665905110;
  solution(165) = -0.1495522317;
  solution(166) = -0.1838028401;
  solution(167) = -0.1954546792;
  solution(168) = -0.1836288060;
  solution(169) = -0.1500743383;
  solution(170) = -0.1669385794;
  solution(171) = -0.1499002842;
  solution(172) = -0.1667645353;
  solution(173) = -0.0411667471;
  solution(174) = -0.0720601654;
  solution(175) = -0.0409926891;
  solution(176) = -0.1029536023;
  solution(177) = -0.0718861173;
  solution(178) = -0.1033017306;
  solution(179) = -0.1266880247;
  solution(180) = -0.1031276566;
  solution(181) = -0.1265139607;
  solution(182) = -0.0415149152;
  solution(183) = -0.0724083136;
  solution(184) = -0.0413408213;
  solution(185) = -0.0722342296;
  solution(186) = -0.1263399164;
  solution(187) = -0.1027795642;
  solution(188) = -0.1261658882;
  solution(189) = -0.1859582178;
  solution(190) = -0.1894319687;
  solution(191) = -0.1965323680;
  solution(192) = -0.1594437570;
  solution(193) = -0.1594436879;
  solution(194) = -0.1761747383;
  solution(195) = -0.1524962754;
  solution(196) = -0.1559700063;
  solution(197) = -0.1692272367;
  solution(198) = -0.1727009777;
  solution(199) = -0.1102697214;
  solution(200) = -0.1102696718;
  solution(201) = -0.1348566702;
  solution(202) = -0.0420842142;
  solution(203) = -0.0420841843;
  solution(204) = -0.0761769187;
  solution(205) = -0.0420843333;
  solution(206) = -0.0420842638;
  solution(207) = -0.0761770180;
  solution(208) = -0.0761769584;
  solution(209) = -0.1033224380;
  solution(210) = -0.1067961490;
  solution(211) = -0.1279093469;
  solution(212) = -0.1102698797;
  solution(213) = -0.1313830679;
  solution(214) = -0.0420845314;
  solution(215) = -0.0420844225;
  solution(216) = -0.0761771961;
  solution(217) = -0.0761770972;
  solution(218) = -0.0351371296;
  solution(219) = -0.0386108207;
  solution(220) = -0.0692297745;
  solution(221) = -0.0727034754;
  solution(222) = -0.1102697907;
  solution(223) = -0.1348568085;
  solution(224) = -0.1348567295;
  solution(225) = -0.0203667548;
  solution(226) = -0.0512601332;
  solution(227) = -0.0309408252;
  solution(228) = -0.0821535303;
  solution(229) = -0.0618342136;
  solution(230) = -0.0486917077;
  solution(231) = -0.0720779618;
  solution(232) = -0.0654226091;
  solution(233) = -0.0954642355;
  solution(234) = -0.0888088733;
  solution(235) = 0.0130950279;
  solution(236) = -0.0177983306;
  solution(237) = -0.0036358536;
  solution(238) = -0.0345292221;
  solution(239) = -0.0462904780;
  solution(240) = -0.0631546592;
  solution(241) = -0.0708773470;
  solution(242) = -0.0800188600;
  solution(243) = -0.0877415381;
  solution(244) = -0.0118335915;
  solution(245) = -0.0234853507;
  solution(246) = -0.0459262164;
  solution(247) = -0.0575779856;
  solution(248) = 0.0218947505;
  solution(249) = 0.0050305893;
  solution(250) = -0.0121978544;
  solution(251) = -0.0290620256;
  solution(252) = 0.0622687055;
  solution(253) = 0.0313753669;
  solution(254) = 0.0376818764;
  solution(255) = 0.0004820098;
  solution(256) = 0.0067885279;
  solution(257) = 0.0686671985;
  solution(258) = 0.0452809842;
  solution(259) = 0.0345746135;
  solution(260) = 0.0111883893;
  solution(261) = 0.1304538543;
  solution(262) = 0.0995605357;
  solution(263) = 0.0963612892;
  solution(264) = 0.0654679606;
  solution(265) = -0.0229042244;
  solution(266) = -0.0241048392;
  solution(267) = -0.0474910833;
  solution(268) = -0.0916706392;
  solution(269) = -0.1291926574;
  solution(270) = -0.1123284366;
  solution(271) = -0.1046057489;
  solution(272) = -0.1408444565;
  solution(273) = -0.1162575381;
  solution(274) = -0.1289260980;
  solution(275) = -0.1055398044;
  solution(276) = -0.1121951569;
  solution(277) = -0.0927276206;
  solution(278) = -0.1395002083;
  solution(279) = -0.1161139047;
  solution(280) = -0.1743063790;
  solution(281) = -0.1626545598;
  solution(282) = -0.1575754079;
  solution(283) = -0.1732286901;
  solution(284) = -0.1563644394;
  solution(285) = -0.1848805192;
  solution(286) = -0.1457903191;
  solution(287) = -0.1459235988;
  solution(288) = -0.1290593680;
  solution(289) = -0.0713640918;
  solution(290) = -0.0404706338;
  solution(291) = -0.0402966551;
  solution(292) = -0.1026055459;
  solution(293) = -0.0717120890;
  solution(294) = -0.0715380805;
  solution(295) = -0.0408186508;
  solution(296) = -0.0406446324;
  solution(297) = -0.1259918798;
  solution(298) = -0.1024315473;
  solution(299) = -0.1029536023;
  solution(300) = -0.0720601654;
  solution(301) = -0.0718861173;
  solution(302) = -0.0411667471;
  solution(303) = -0.0409926891;
  solution(304) = -0.1033017306;
  solution(305) = -0.0724083136;
  solution(306) = -0.0722342296;
  solution(307) = -0.0415149152;
  solution(308) = -0.0413408213;
  solution(309) = -0.1500743383;
  solution(310) = -0.1266880247;
  solution(311) = -0.1265139607;
  solution(312) = -0.1031276566;
  solution(313) = -0.1665905110;
  solution(314) = -0.1497262499;
  solution(315) = -0.1495522317;
  solution(316) = -0.1838028401;
  solution(317) = -0.1669385794;
  solution(318) = -0.1667645353;
  solution(319) = -0.1499002842;
  solution(320) = -0.1954546792;
  solution(321) = -0.1836288060;
  solution(322) = -0.1263399164;
  solution(323) = -0.1261658882;
  solution(324) = -0.1027795642;
  solution(325) = -0.0821535303;
  solution(326) = -0.0512601332;
  solution(327) = -0.0618342136;
  solution(328) = -0.0203667548;
  solution(329) = -0.0309408252;
  solution(330) = -0.0486917077;
  solution(331) = -0.0177983306;
  solution(332) = -0.0345292221;
  solution(333) = 0.0130950279;
  solution(334) = -0.0036358536;
  solution(335) = -0.0954642355;
  solution(336) = -0.0720779618;
  solution(337) = -0.0888088733;
  solution(338) = -0.0654226091;
  solution(339) = 0.0004820098;
  solution(340) = 0.0313753669;
  solution(341) = 0.0067885279;
  solution(342) = 0.0622687055;
  solution(343) = 0.0376818764;
  solution(344) = 0.0686671985;
  solution(345) = 0.0995605357;
  solution(346) = 0.0654679606;
  solution(347) = 0.1304538543;
  solution(348) = 0.0963612892;
  solution(349) = 0.0218947505;
  solution(350) = 0.0452809842;
  solution(351) = 0.0111883893;
  solution(352) = 0.0345746135;
  solution(353) = -0.0800188600;
  solution(354) = -0.0631546592;
  solution(355) = -0.0877415381;
  solution(356) = -0.0462904780;
  solution(357) = -0.0708773470;
  solution(358) = -0.0118335915;
  solution(359) = 0.0050305893;
  solution(360) = -0.0290620256;
  solution(361) = -0.0121978544;
  solution(362) = -0.0351371296;
  solution(363) = -0.0234853507;
  solution(364) = -0.0575779856;
  solution(365) = -0.0459262164;
  solution(366) = -0.0229042244;
  solution(367) = -0.0474910833;
  solution(368) = -0.0241048392;
  solution(369) = -0.1894319687;
  solution(370) = -0.1859582178;
  solution(371) = -0.1965323680;
  solution(372) = -0.1594437570;
  solution(373) = -0.1559700063;
  solution(374) = -0.1727009777;
  solution(375) = -0.1524962754;
  solution(376) = -0.1692272367;
  solution(377) = -0.1594436879;
  solution(378) = -0.1761747383;
  solution(379) = -0.1102698797;
  solution(380) = -0.1067961490;
  solution(381) = -0.1313830679;
  solution(382) = -0.1033224380;
  solution(383) = -0.1279093469;
  solution(384) = -0.0420845314;
  solution(385) = -0.0386108207;
  solution(386) = -0.0727034754;
  solution(387) = -0.0692297745;
  solution(388) = -0.0420843333;
  solution(389) = -0.0420844225;
  solution(390) = -0.0761770972;
  solution(391) = -0.0761771961;
  solution(392) = -0.1102696718;
  solution(393) = -0.1102697214;
  solution(394) = -0.1348566702;
  solution(395) = -0.0420842142;
  solution(396) = -0.0420842638;
  solution(397) = -0.0761769584;
  solution(398) = -0.0761770180;
  solution(399) = -0.0420841843;
  solution(400) = -0.0761769187;
  solution(401) = -0.1102697907;
  solution(402) = -0.1348567295;
  solution(403) = -0.1348568085;
  solution(404) = -0.0916706392;
  solution(405) = -0.1291926574;
  solution(406) = -0.1408444565;
  solution(407) = -0.1162575381;
  solution(408) = -0.1123284366;
  solution(409) = -0.1046057489;
  solution(410) = -0.1626545598;
  solution(411) = -0.1743063790;
  solution(412) = -0.1575754079;
  solution(413) = -0.1848805192;
  solution(414) = -0.1563644394;
  solution(415) = -0.1732286901;
  solution(416) = -0.1055398044;
  solution(417) = -0.1289260980;
  solution(418) = -0.1121951569;
  solution(419) = -0.1161139047;
  solution(420) = -0.1395002083;
  solution(421) = -0.0927276206;
  solution(422) = -0.1457903191;
  solution(423) = -0.1290593680;
  solution(424) = -0.1459235988;
  solution(425) = 0.0730639635;
  solution(426) = 0.0730639720;
  solution(427) = 0.1041312494;
  solution(428) = 0.0730639635;
  solution(429) = 0.1041312494;
  solution(430) = 0.1351985169;
  solution(431) = 0.0109293702;
  solution(432) = 0.0184365031;
  solution(433) = 0.0419966761;
  solution(434) = 0.0259436179;
  solution(435) = 0.0495038004;
  solution(436) = 0.0259436179;
  solution(437) = 0.0184365031;
  solution(438) = 0.0495038004;
  solution(439) = 0.0419966761;
  solution(440) = -0.0361910153;
  solution(441) = -0.0221618094;
  solution(442) = -0.0051236894;
  solution(443) = -0.0081326219;
  solution(444) = 0.0089055078;
  solution(445) = -0.0702672950;
  solution(446) = -0.0510256673;
  solution(447) = -0.0391999491;
  solution(448) = -0.0317840580;
  solution(449) = -0.0199583301;
  solution(450) = -0.0833114407;
  solution(451) = -0.0767893583;
  solution(452) = -0.0597512182;
  solution(453) = -0.0532291453;
  solution(454) = -0.0081326219;
  solution(455) = -0.0221618094;
  solution(456) = 0.0089055078;
  solution(457) = -0.0361910153;
  solution(458) = -0.0051236894;
  solution(459) = -0.0702672950;
  solution(460) = -0.0767893583;
  solution(461) = -0.0532291453;
  solution(462) = -0.0597512182;
  solution(463) = -0.0317840580;
  solution(464) = -0.0510256673;
  solution(465) = -0.0199583301;
  solution(466) = -0.0391999491;
  solution(467) = -0.0361910058;
  solution(468) = -0.0126308128;
  solution(469) = -0.0126308128;
  solution(470) = -0.0939187710;
  solution(471) = -0.0664990555;
  solution(472) = -0.0628514052;
  solution(473) = -0.0390793584;
  solution(474) = -0.0354316983;
  solution(475) = -0.1012141113;
  solution(476) = -0.0703205760;
  solution(477) = -0.0701467255;
  solution(478) = -0.0394270591;
  solution(479) = -0.0392531989;
  solution(480) = -0.1410392763;
  solution(481) = -0.1211266842;
  solution(482) = -0.1174790139;
  solution(483) = -0.0975664313;
  solution(484) = -0.1015618519;
  solution(485) = -0.0706683357;
  solution(486) = -0.0704944462;
  solution(487) = -0.0397748380;
  solution(488) = -0.0396009387;
  solution(489) = -0.0708422450;
  solution(490) = -0.0399487572;
  solution(491) = -0.1251221347;
  solution(492) = -0.1017357514;
  solution(493) = -0.1751156759;
  solution(494) = -0.1617251565;
  solution(495) = -0.1580774663;
  solution(496) = -0.1483346565;
  solution(497) = -0.1446869566;
  solution(498) = -0.1653728665;
  solution(499) = -0.1485085369;
  solution(500) = -0.1869414640;
  solution(501) = -0.1787633761;
  solution(502) = -0.1249482446;
  solution(503) = -0.1247743741;
  solution(504) = -0.1013879717;
  solution(505) = -0.0390793584;
  solution(506) = -0.0664990555;
  solution(507) = -0.0354316983;
  solution(508) = -0.0939187710;
  solution(509) = -0.0628514052;
  solution(510) = -0.1012141113;
  solution(511) = -0.1211266842;
  solution(512) = -0.0975664313;
  solution(513) = -0.1410392763;
  solution(514) = -0.1174790139;
  solution(515) = -0.0394270591;
  solution(516) = -0.0703205760;
  solution(517) = -0.0392531989;
  solution(518) = -0.0701467255;
  solution(519) = -0.1483346565;
  solution(520) = -0.1617251565;
  solution(521) = -0.1446869566;
  solution(522) = -0.1751156759;
  solution(523) = -0.1580774663;
  solution(524) = -0.1787633761;
  solution(525) = -0.1869414640;
  solution(526) = -0.1485085369;
  solution(527) = -0.1653728665;
  solution(528) = -0.0397748380;
  solution(529) = -0.0706683357;
  solution(530) = -0.0396009387;
  solution(531) = -0.1015618519;
  solution(532) = -0.0704944462;
  solution(533) = -0.1017357514;
  solution(534) = -0.1251221347;
  solution(535) = -0.0399487572;
  solution(536) = -0.0708422450;
  solution(537) = -0.1249482446;
  solution(538) = -0.1013879717;
  solution(539) = -0.1247743741;
  solution(540) = -0.1751156662;
  solution(541) = -0.1514641200;
  solution(542) = -0.1462516885;
  solution(543) = -0.1632898881;
  solution(544) = -0.1462516885;
  solution(545) = -0.1632898881;
  solution(546) = -0.1173877604;
  solution(547) = -0.1056532561;
  solution(548) = -0.1292135085;
  solution(549) = -0.0820930231;
  solution(550) = -0.1003495907;
  solution(551) = -0.0938275179;
  solution(552) = -0.1056532561;
  solution(553) = -0.1173877604;
  solution(554) = -0.1292135085;
  solution(555) = -0.0938275179;
  solution(556) = -0.1003495907;
  solution(557) = -0.0820930231;
  solution(558) = -0.1173877507;
  solution(559) = -0.1344259303;
  solution(560) = -0.1344259303;

  Eigen::VectorXd DG_solution(dofs_DG);

  //convert to DG solution
  solver.c0_converter.init(solver.grid, DG_solution.size());
  solver.c0_converter.convert_coefficients_toDG(solution, DG_solution);

  // write solution into leaf cells
  solver.restore_MA(DG_solution);

  //plot solution to file *Inter.vtu
  solver.plotter.write_numericalsolution_VTK(0, true);

#else
  cerr << "Unknown equation! Program aborted!" << endl;
#endif


  return 0;

}
