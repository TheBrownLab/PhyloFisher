/*
 * Divvier.cpp
 *
 *  Created on: 24 Oct 2017
 *      Author: simon
 *  ---
 *  A program for divvying or partial filtering of alignments
 *
 *  Usage:
 *   ./divvier [OPTIONS] file threshold
 *
 *  Options choosing clustering approach:
 *  -divvy 		  : do standard divvying (DEFAULT)
 *  -partial 	  : Do only partial filtering by checking whether each character is joined with the rest
 *  (-partialall) : do the partial approach, but output the characters from the remaining divvied columns as single characters. Needed for MSA accuracy calcs
 *
 * NOTE: Experimentally derived thresholds are: 0.857 for partial filtering (FDR = 0.1) and 0.801 for divvying (FDR = 0.01)
 *
 *
 *  Options for defining the heuristics
 *  -approx X     : Minimum number of comparisons required for testing each partition in the divvied MSA (DEFAULT: X = 10; if X < 0 then do all)
 *  -HMMapprox    : Do the pairHMM bounding approximation (DEFAULT)
 *  -HMMexact     : Do the full pairHMM and ignore bounding
 */

// Imported stuff from Zorro and bionj
extern "C" {
#define EXTERN extern
#include "hmm.h"
#include "utils.h"
#include "matrices.h"
#undef EXTERN

}

#include "Divvier.h"
#include "Cluster.h"
#include "Tree.h"
#include "Sequence.h"
#include <fstream>

string BigTree = "((EOG52FR0S|AGAMB_3.6|Bombylius_major|s3305_L_37822_0-179:0.10258365074488022539,((EOG52FR0S|AGAMB_3.6|Trichocera_fuscata|C1011097-176:0.37919204223948976828,EOG52FR0S|AGAMB_3.6|Bibio_marci|s6005_L_37988_0-160:0.21668222714464321910):0.07156413693204413673,((EOG52FR0S|NVITR_1.2|Apachyus_chartaceus|s5108_L_54582_0-173:0.27997964631343663644,EOG52FR0S|TCAST_3.0reorderedheaderfields|Forficula_auricularia|C894882-184:0.35947113088008203485):0.40581948723979877069,((((EOG52FR0S|ZNEVA_2.1|Platycentropus_radiatus|s10078_L_106168_0-182:0.28712730117582346834,EOG52FR0S|BMORI_2.0|Rhyacophila_fasciata|C295776-176:0.20879844376069728318):0.08377470818577073541,(((((EOG52FR0S|BMORI_2.0|BMORI|BGIBMGA013243_TA-192:0.20796019754844122240,(EOG52FR0S|BMORI_2.0|Parides_arcas|s11100_L_89186_0-191:0.15948176694983298707,(EOG52FR0S|BMORI_2.0|Polyommatus_icarus|C448597-191:0.20728969418089865373,EOG52FR0S|BMORI_2.0|Zygaena_fausta|C408340-190:0.31914552697676795701):0.07335145172361670629):0.05053893709918481220):0.05850196533057449438,EOG52FR0S|BMORI_2.0|Yponomeuta_evonymella|C659571-189:0.24661963489412574990):0.04379433474415747596,(EOG52FR0S|BMORI_2.0|Nemophora_deegerella|s9473_L_47391_0-192:0.34388825491603491891,(EOG52FR0S|BMORI_2.0|Eriocrania_subpurpurella|C657810-162:0.28016145828226290959,EOG52FR0S|BMORI_2.0|Triodia_sylvina|s19689_L_226945_0-114:0.25180531327347416282):0.05477461074677487246):0.05850214087699501936):0.27893477336102578956,EOG52FR0S|TCAST_3.0reorderedheaderfields|Hydroptilidae_sp|C371109-171:0.20954375979877706837):0.08227257382311191358,(EOG52FR0S|AGAMB_3.6|Philopotamus_ludificatus|contig11196-112:0.18670733170122680300,EOG52FR0S|AGAMB_3.6|Cheumatopsyche_sp|C346787-143:0.24975404489548827525):0.10840344418213380961):0.03058718575931992922):0.20521928397267710786,(((EOG52FR0S|ZNEVA_2.1|Ceratophyllus_gallinae|C304111-181:0.15468533185857691326,(EOG52FR0S|ZNEVA_2.1|Ctenocephalides_felis|C330557-186:0.06671445089985969523,EOG52FR0S|ZNEVA_2.1|Archaeopsylla_erinacei|contig19363-160:0.07554576633173239186):0.09865183937584097451):0.15539737916237764126,(EOG52FR0S|AGAMB_3.6|Bittacus_pilicornis|C194995-166:0.13838638337865413752,EOG52FR0S|TCAST_3.0reorderedheaderfields|Panorpa_vulgaris|C326395-171:0.19468361452809271328):0.08887788353119491225):0.06652652542145692793,(EOG52FR0S|DMELA_5.40|Nannochorista_sp|contig01540-115:0.31122981729268078821,EOG52FR0S|ZNEVA_2.1|Boreus_hyemalis|s9185_L_107660_0-179:0.20765950329740226477):0.05863533568266944551):0.03902888353479667255):0.03538387372803234593,(((((EOG52FR0S|ZNEVA_2.1|Zorotypus_caudelli|C2490937-177:0.34715128115262866570,(EOG52FR0S|ZNEVA_2.1|Grylloblatta_bifratrilecta|s13116_L_117268_0-182:0.04666403811736503232,EOG52FR0S|ZNEVA_2.1|Galloisiana_yuasai|C897547-182:0.08611999898518697683):0.26871467228598600041):0.03144048514876583017,((((EOG52FR0S|ZNEVA_2.1|Ceuthophilus_sp|C1389194-183:0.24510371778815792654,(EOG52FR0S|ZNEVA_2.1|Acanthocasuarina_muellerianae|C679096-177:0.78333252838140332575,EOG52FR0S|ZNEVA_2.1|Gryllotalpa_sp|s2139_L_10561_0-181:0.30716835037757195259):0.06356106894554112985):0.12891568640354211794,(EOG52FR0S|ZNEVA_2.1|Nilaparvata_lugens|C572437-185:0.30770644845472772122,((EOG52FR0S|ZNEVA_2.1|Timema_cristinae|s15504_L_191927_0-183:0.31245666328246152199,((EOG52FR0S|ZNEVA_2.1|Campodea_augens|s12453_L_61651_1-175:0.64320360852185343159,EOG52FR0S|ZNEVA_2.1|Occasjapyx_japonicus|C446479-170:0.63303573982941852005):0.40917601625516664132,(EOG52FR0S|AGAMB_3.6|Velia_caprai|s10579_L_125610_0-173:0.66117482127068250009,((EOG52FR0S|ZNEVA_2.1|Peruphasma_schultei|s5798_L_35856_0-182:0.11847229589514635117,EOG52FR0S|ZNEVA_2.1|Aretaon_asperrimus|s17907_L_301633_0-182:0.07253143523751401367):0.19903205412866231683,EOG52FR0S|ZNEVA_2.1|Ranatra_linearis|s1330_L_4194_1-169:0.55354515933741088585):0.03262464634031971705):0.04083863261130472183):0.17194594463365361903):0.13632695045330456285,(EOG52FR0S|ZNEVA_2.1|Haploembia_palaui|C892465-181:0.03062609439495739272,EOG52FR0S|ZNEVA_2.1|Aposthonia_japonica|C1112469-182:0.08760650591975023549):0.20538604134669566359):0.06202511515297982891):0.03988691757109619901):0.06753599270869413418,((EOG52FR0S|ZNEVA_2.1|Tetrix_subulata|C601439-182:0.37131693579780711278,(EOG52FR0S|ZNEVA_2.1|Prosarthria_teretrirostris|C772607-153:0.27003519586186774948,EOG52FR0S|ZNEVA_2.1|Stenobothrus_lineatus|C2019176-167:0.31558548049701179439):0.08907778373876351630):0.07790452599305680570,(EOG52FR0S|ZNEVA_2.1|Perla_marginata|C644184-178:0.20904642311162455193,(EOG52FR0S|ZNEVA_2.1|Cosmioperla_kuna|C572547-181:0.28459243504544345926,EOG52FR0S|ZNEVA_2.1|Leuctra_sp|C439181-170:0.23534844734160628721):0.03611738912271651725):0.05315439554639348613):0.00628857273714119279):0.03777628665101530336,((EOG52FR0S|ZNEVA_2.1|Metallyticus_splendidus|s21420_L_361771_0-182:0.08139199083812696800,(EOG52FR0S|ZNEVA_2.1|Empusa_pennata|C1257614-182:0.05831369988802588555,EOG52FR0S|ZNEVA_2.1|Mantis_religiosa|s16763_L_273691_0-182:0.03767170927567057431):0.06640196255228469902):0.22577155249681363225,(EOG52FR0S|ZNEVA_2.1|Blaberus_atropus|s21873_L_207274_0-184:0.14385889001104165685,(EOG52FR0S|ZNEVA_2.1|Cryptocercus_sp|C1038622-184:0.06910321542211660117,((EOG52FR0S|ZNEVA_2.1|Periplaneta_americana|C1181093-183:0.19336899542727936652,EOG52FR0S|ZNEVA_2.1|Mastotermes_darwiniensis|C2113619-184:0.05240389976927931764):0.02280219262590709936,(EOG52FR0S|ZNEVA_2.1|ZNEVA|Znev_01260-184:0.07448774569891650210,EOG52FR0S|ZNEVA_2.1|Prorhinotermes_simplex|s11217_L_143197_0-184:0.06831791999531773574):0.00531771294730884643):0.02772591341670576498):0.13298325015205483823):0.05305481403293376558):0.08813889622810555657):0.01032094283925216564):0.04918804217134938039,((EOG52FR0S|ZNEVA_2.1|Gynaikothrips_ficorum|C875677-171:0.34817770995555730185,(EOG52FR0S|ZNEVA_2.1|Thrips_palmi|C276765-175:0.15906098279260383332,EOG52FR0S|TCAST_3.0reorderedheaderfields|Frankliniella_cephalica|C380329-158:0.21388659633192913523):0.16048147894155415094):0.07653538013830162023,((((EOG52FR0S|ZNEVA_2.1|Cypridininae_sp|DMPC15491392-165:0.52790443634768879910,(EOG52FR0S|PHUMA_1.2|Ectopsocus_briggsi|s15732_L_91826_0-174:0.46448893123991685794,(EOG52FR0S|PHUMA_1.2|Liposcelis_bostrychophila|C576092-177:0.30197120987225656297,EOG52FR0S|PHUMA_1.2|PHUMA|PHUM000230_RA-190:0.12517489564931963408):0.09559427078257826116):0.13402161218995289893):0.06452847926333950268,(EOG52FR0S|ZNEVA_2.1|Thermobia_domestica|C1632884-171:0.28087239716863043881,(EOG52FR0S|AECHI_3.8|Lepeophtheirus_salmonis|347703492-168:2.69374416389131710048,EOG52FR0S|DPULE_2010beta3|DPULE|hxNCBI_GNO_158323-183:0.39413559219271232514):0.26879970332748343020):0.03490067545383494785):0.03655927293058804789,((EOG52FR0S|PHUMA_1.2|Planococcus_citri|C485257-176:0.57727675594194338693,(EOG52FR0S|ZNEVA_2.1|Pogonognathellus_lon_fla|s10192_L_40563_0-168:0.35212892139224716281,((EOG52FR0S|AGAMB_3.6|Tetrodontophora_bielanensis|s4398_L_11904_0-170:0.31253292586648845353,EOG52FR0S|ZNEVA_2.1|Anurida_maritima|s7480_L_41346_1-179:0.53534165142180289987):0.09780428551627755318,(EOG52FR0S|PHUMA_1.2|Folsomia_candida|s13407_L_107771_0-183:0.32730306719440233065,EOG52FR0S|ZNEVA_2.1|Sminthurus_vir_nig|C438893-180:0.44815131648704759071):0.09981773077066337374):0.04633260513791208346):0.14420707809006821920):0.07789073415241687393,(((EOG52FR0S|ISCAP_1.1|ISCAP|ISCW010767_RA-187:0.89808705782006736928,EOG52FR0S|ZNEVA_2.1|Glomeris_pustulata|DMPC15855832-177:0.26687653412383394169):0.13569990008012094984,EOG52FR0S|ZNEVA_2.1|Symphylella_vulgaris|DMPC15369782-167:0.18976574371435883659):0.32321515569392572642,(EOG52FR0S|ZNEVA_2.1|Meinertellus_cundinamarcensis|C1499868-174:0.24727594137049682677,EOG52FR0S|ZNEVA_2.1|Machilis_hrabei|C992353-174:0.27952995969771998741):0.05531855326546314400):0.09807100926610833047):0.00000284652651879708):0.06328495355362576125,(((((EOG52FR0S|DPULE_2010beta3|Ephemera_danica|C517837-183:0.26107317656135214934,EOG52FR0S|DPULE_2010beta3|Eurylophella_sp|s1786_L_4814_0-127:0.22435150151057325907):0.26812301515613495839,EOG52FR0S|ZNEVA_2.1|Isonychia_bicolor|s3662_L_10128_0-182:0.21626403668416263604):0.09224864577838710888,(EOG52FR0S|ZNEVA_2.1|Acerentomon_sp|C729640-170:0.42360344963029672449,(EOG52FR0S|TCAST_3.0reorderedheaderfields|Litopenaeus_vannamei|DMPC7768585-59:0.92986603079020624385,EOG52FR0S|ZNEVA_2.1|Baetis_pumilus|C418006-172:0.24953913442093139663):0.06252057145089949530):0.09415038571919742694):0.12507299422017811863,((EOG52FR0S|ZNEVA_2.1|Speleonectes_tulumensis|DMPC15231531-108:0.45393692433976629008,EOG52FR0S|ZNEVA_2.1|Tricholepidion_gertschi|C1162185-171:0.16794054281747641810):0.05393389030309596321,EOG52FR0S|ZNEVA_2.1|Atelura_formicaria|C1786565-181:0.29696880472858544486):0.01849951439992951474):0.07567091885262067219,((EOG52FR0S|ZNEVA_2.1|Cordulegaster_boltoni|C513443-179:0.08046521448541195387,EOG52FR0S|ZNEVA_2.1|Epiophlebia_superstes|s9066_L_135880_0-171:0.05969989235595130755):0.08201709350185917846,EOG52FR0S|ZNEVA_2.1|Calopteryx_splendens|s6555_L_128950_0-175:0.18384783375122812354):0.19182395594675868966):0.01807997043484640964):0.01795310546278424194):0.06409461969434167294):0.05483711695084447085,(((EOG52FR0S|ZNEVA_2.1|Okanagana_villosa|C1047067-175:0.26389540938051336827,EOG52FR0S|ZNEVA_2.1|Tanzaniophasma_sp|C698355-178:0.34635230878060485615):0.06844350128053054705,(EOG52FR0S|TCAST_3.0reorderedheaderfields|Acanthosoma_haemorrhoidale|s6806_L_55372_0-167:0.41505233401759783485,EOG52FR0S|AGAMB_3.6|Xenophysella_greensladeae|s20205_L_145124_0-178:0.34376572523914900037):0.09259141827950416459):0.10450190963043513859,(((((EOG52FR0S|TCAST_3.0reorderedheaderfields|Xanthostigma_xanthostigma|C566271-164:0.19449936655034971711,EOG52FR0S|PHUMA_1.2|Inocellia_crassicornis|s10337_L_80389_0-154:0.33767557295784456084):0.07652478674593690688,((EOG52FR0S|TCAST_3.0reorderedheaderfields|Osmylus_fulvicephalus|C325502-162:0.32724880512702875235,(EOG52FR0S|TCAST_3.0reorderedheaderfields|Euroleon_nostras|C417689-154:0.22722803678459393972,EOG52FR0S|TCAST_3.0reorderedheaderfields|Dichochrysa_prasina|C865900-157:0.17965260506147759378):0.09226201333368116986):0.04968307450628078881,(EOG52FR0S|AECHI_3.8|Conwentzia_psociformis|s5518_L_8149_2-164:0.90973411417670912993,(EOG52FR0S|TCAST_3.0reorderedheaderfields|Corydalus_cornutus|C168679-171:0.36817967172391224961,EOG52FR0S|TCAST_3.0reorderedheaderfields|Sialis_lutaria|contig05237-171:0.31588037603904894901):0.08098061631711647723):0.02889324885925097230):0.21217546867391573473):0.42075910593195908760,((EOG52FR0S|TCAST_3.0reorderedheaderfields|Priacma_serrata|contig01572-167:0.48448965335818028333,(EOG52FR0S|TCAST_3.0reorderedheaderfields|Aleochara_curtula|C207305-170:0.31950480379657325569,(EOG52FR0S|TCAST_3.0reorderedheaderfields|Dendroctonus_ponderosae|DMPC14949187-172:0.25343465759882993771,(EOG52FR0S|TCAST_3.0reorderedheaderfields|TCAST|TC015202-338:0.39106414960961038974,EOG52FR0S|TCAST_3.0reorderedheaderfields|Meloe_violaceus|s3323_L_15628_0-171:0.32192844068618448050):0.11815592516250769672):0.06869563842199601089):0.02524545323405661193):0.06469464556989516779,(EOG52FR0S|TCAST_3.0reorderedheaderfields|Carabus_granulatus|contig20817-171:0.20936975200398974528,(EOG52FR0S|TCAST_3.0reorderedheaderfields|Lepicerus_sp|s31954_L_190160_0-169:0.30174290416189342157,(EOG52FR0S|TCAST_3.0reorderedheaderfields|Gyrinus_marinus|C426169-130:0.31698713369864911504,(EOG52FR0S|TCAST_3.0reorderedheaderfields|Stylops_melittae|s815_L_3558_4-173:1.09132541448030906395,EOG52FR0S|BMORI_2.0|Mengenilla|contig00812-50:0.24860954514891062117):0.36817374316011364233):0.09639351752230840287):0.12968634453582394972):0.04153894268180059052):0.09005386623017200276):0.07505777768422575158,EOG52FR0S|TCAST_3.0reorderedheaderfields|Cercopis_vulnerata|s2282_L_9869_0-164:0.44816506587907789516):0.04451910789939503288,(EOG52FR0S|TCAST_3.0reorderedheaderfields|Bemisia_tabaci|gi321197821-156:0.18512744296802555177,EOG52FR0S|TCAST_3.0reorderedheaderfields|Trialeurodes_vaporariorum|C822518-167:0.16601394260731558439):0.28381702924840379598):0.06874419825445349241):0.01224687049888137542):0.03621926839709280199,(EOG52FR0S|AECHI_3.8|Cotesia_vestalis|C300420-176:0.35530433915661974176,((EOG52FR0S|NVITR_1.2|NVITR|NV18101_RA-206:0.25865233882925381392,((EOG52FR0S|APISU_v2|Essigella_californica|s21039_L_220779_0-183:0.21768785384887601175,EOG52FR0S|APISU_v2|APISU|ACYPI006767_RA-183:0.00000284652651879708):0.97560016268775173742,EOG52FR0S|NVITR_1.2|Leptopilina_clavipes|C295537-168:0.21594239260245212675):0.05378157473385462156):0.07305057817828163047,(EOG52FR0S|AECHI_3.8|Chrysis_viridula|s10150_L_87082_0-168:0.18506806268268455318,(EOG52FR0S|AECHI_3.8|Tenthredo_koehleri|C329071-193:0.29226817370556612552,((EOG52FR0S|AECHI_3.8|Exoneura_robusta|gi323881566-124:0.16346549139818158936,EOG52FR0S|AMELL_prerelease2|AMELI|GB13139_RA-177:0.25344702471305829983):0.09054987763675816093,(EOG52FR0S|AECHI_3.8|Harp_salt|Hsal_15250-127:0.08760556090348606273,EOG52FR0S|AECHI_3.8|AECHI|Aech_11646-193:0.17994232010852431736):0.08338187287694917571):0.05369313507965074034):0.02027463668964252577):0.10887042098547439206):0.01448637001041356773):0.20399113308697866542):0.03794794241833907011):0.07883567945606387295):0.06444902158795283442):0.07803962599227282082):0.16751648001610505712,(((EOG52FR0S|AGAMB_3.6|Notostira_elongata|C433581-166:0.71271863719522066116,(EOG52FR0S|AGAMB_3.6|AGAMB|AGAP003199_RA-204:0.15983289722980076331,EOG52FR0S|AGAMB_3.6|Aede_aegy|AAEL000163-RA-202:0.16344503318490347099):0.02628487693242425455):0.26383259074490533758,EOG52FR0S|DMELA_5.40|DMELA|3674_CG31229PA-195:0.33554437029363948231):0.15838475752009065212,(EOG52FR0S|DMELA_5.40|Rhagoletis_pomonella|gi257100020-155:0.34539655929732071549,(EOG52FR0S|DMELA_5.40|Sarcophaga_crassipalpis|gi296337228-179:0.08465742730590501697,EOG52FR0S|DMELA_5.40|Triarthria_setipennis|C437403-177:0.14214045150649035065):0.07479398006429255341):0.04071023310863794431):0.07834467316795921954,EOG52FR0S|DMELA_5.40|Lipara_lucens|C642472-194:0.25694883168138404894);";

vector <CPostP> allPP;		// Global holding the posteriors
vector <string> names;
vector <string> in_seq;
void CreateZorro(); // Must be called after names and in_seq are initialised
string MakeTreeNJDist(std::vector <std::string> Names, std::vector <std::string> Seqs);	// Function that actually builds the tree

const double divvyThreshold = 0.801;		// 1% FDR
const double partialThreshold = 0.774;		// Needs confirming, but currently 10% FDR due to sensitivity error

// Structure defining options
struct SOptions {
	// File information
	string suffixPP = ".PP";
	string suffixDivvy = ".divvy";
	// Divvying information
	int doFullDivvy = 1; 				// 1 = standard divvying; 0 = do partial filtering; -1 do partial filtering with messy output
	bool HMMapproximation = true;		// Whether to perform the pairHMM approximation
	bool acceptNoInfo = false;			// Whether comparisons between sets when divvying are accepted if there's no information
	int approxNumber = 10;
	double threshold = divvyThreshold;
	bool forceValidate = false;
	// Output info
	int minDivvyChar = 2;
	string divvy_char = "*";
} options;

int main(int argc, char *argv[]) {
	int colCount = 0;
	string fileIn;
	vector <CSequence> *seqData = NULL;

	cout << "\n========================================================================";
	cout << "\n    Divvier ("<<VERSION_NUMBER<<"): a program for MSA processing by Simon Whelan";
	cout << "\n========================================================================";

	bool showHelp = false;
	bool doneThreshold = false;
	// Process options
	if(argc >= 2) {
		if(strcmp(argv[1], "-h") != 0) {
			// Get the file name and read it
			fileIn = argv[argc - 1];
			seqData = ReadSequences(argv[argc - 1]);
			// Get the options
			for(int i = 1 ; i < argc - 1; i++) {
				if(argv[i][0] != '-') { continue; }
				else if(strcmp(argv[i],"-divvy") == 0) { options.suffixDivvy = ".divvy"; options.doFullDivvy = 1; if(!doneThreshold) { options.threshold = divvyThreshold; } }
				else if(strcmp(argv[i], "-partialall") == 0) { options.suffixDivvy = ".partialall"; options.doFullDivvy = -1; options.threshold = partialThreshold; }
				else if(strcmp(argv[i], "-partial") == 0) {
					options.suffixDivvy = ".partial";
					options.doFullDivvy = 0;
					options.divvy_char = "-";
					if(options.minDivvyChar < 2) { options.minDivvyChar = 2; }
					if(!doneThreshold) { options.threshold = 0.857; }
				}
				else if(strcmp(argv[i], "-HMMapprox") == 0) { options.HMMapproximation = true; }
				else if(strcmp(argv[i], "-HMMexact") == 0) { options.HMMapproximation = false; }
				else if(strcmp(argv[i], "-approx") == 0) {
					if(!CheckNext(++i,argc,argv)) { cout << "\nError when setting X in -approx\n\n"; exit(-1); }
					options.approxNumber = atoi(argv[i]);
					if(!InRange(options.approxNumber,1,100000)) { cout << "\n-approx options requires X to be a positive number\n\n"; exit(-1); }
				}
				else if(strcmp(argv[i], "-thresh") == 0 ) { doneThreshold = true; options.threshold = atof(argv[i+1]); }
				else if(strcmp(argv[i], "-checksplits") == 0) { options.forceValidate = true; }
				else if(strcmp(argv[i], "-mincol") == 0) {
					if(!CheckNext(++i,argc,argv)) { cout << "\nError when setting X in -mincol\n\n"; exit(-1); }
					options.minDivvyChar = atoi(argv[i]);
					if(!InRange(options.minDivvyChar,1,100000)) { cout << "\n-mincol options requires X to be a positive number\n\n"; exit(-1); }
				}
				else if(strcmp(argv[i], "-divvygap") == 0) { options.divvy_char = "-"; }
			}
		} else { showHelp = true; }
	} else { showHelp = true; }
	if(showHelp) {
		cout << "\n\n./divvier [options] input_file ";
		cout << "\n\nClustering options:";
		cout << "\n\t-divvy       : do standard divvying (DEFAULT)";
		cout << "\n\t-partial     : do partial filtering by testing removal of individual characters";
		cout << "\n\t-thresh X    : set the threshold for divvying to X (DEFAULT divvying = " << divvyThreshold << "; partial = " << partialThreshold << ")";
		cout << "\n\nApproximation options: ";
		cout << "\n\t-approx X    : minimum number of characters tested in a split during divvying (DEFAULT X = " << options.approxNumber << ")";
		cout << "\n\t-checksplits : go through sequence and ensure there's a pair for every split. Can be slow";
		cout << "\n\t-HMMapprox   : Do the pairHMM bounding approximation (DEFAULT)";
		cout << "\n\t-HMMexact    : Do the full pairHMM and ignore bounding";
		cout << "\n\nOutput options: ";
		cout << "\n\t-mincol X    : Minimum number of characters in a column to output when divvying/filtering (DEFAULT X = " << options.minDivvyChar << ")";
		cout << "\n\t-divvygap    : Output a gap instead of the static * character so divvied MSAs can be used in phylogeny program";
		cout << "\n\n";
		exit(-1);
	}

	// Apply options
	CCluster::Instance()->SetOptions(options.acceptNoInfo, options.approxNumber, (options.doFullDivvy != 1)?true:false, options.forceValidate);
	DO_HMM_APPROX = options.HMMapproximation;
	if(!InRange(options.threshold,0.0,1.0)) { cout << "\nError: threshold must be in range (0,1)\n"; exit(-1); }
	assert(seqData != NULL);

	// Read in the sequence using CSequence's method and transfer it to the variables that Zorro has in global scope here
	cout << "\nInitialising sequences and pairHMM" << flush;
	for(int i = 0; i < seqData->size() ; i++) {
		names.push_back(seqData->at(i).Name());
		in_seq.push_back(seqData->at(i).Seq());
	}
	CreateZorro();
	cout << " ... done" << flush;
	cout << "\nThere are " << seqData->size() << " sequences in an alignment of length " << alen << flush;

	// Initialise the Clustering program
	CCluster::Instance()->AddNames(names);
	cout << "\nMaking guide tree for splitting: " << flush;
	CCluster::Instance()->AddTree(MakeTree(names,in_seq),in_seq);
	cout << " ... done" << flush;

	// Get the posteriors, either through computation or through HMM calcs
	cout << "\n\nGetting posterior probabilities from pairHMM" << flush;
	if(options.HMMapproximation) { cout << " (approx)" << flush; } else { cout << " (exact)" << flush; }
	GetPosteriors(fileIn + options.suffixPP);

	cout << "\nPerforming ";
	switch(options.doFullDivvy) {
	case -1: cout << "partial filtering  with full output"; break;
	case 0:  cout << "partial filtering"; break;
	case 1:  cout << "full divvying"; break;
	default:
		cout << "\nUnknown options for divvying..."; exit(-1);
	}
	cout << " along alignment with threshold " << options.threshold;
	cout << "\nEach column required to have " << options.minDivvyChar << " characters" << endl << flush;

	vector <stringstream> out_seq(Nseq);
	vector <double> PP(Nseq*Nseq,1);

	// Do the divvying and output
	for(int i = 0; i < alen; i++) {
		ProgressSpinner(i+1,alen,"\t");
		// Get the appropriate PPs
		for(auto & x : PP) { x = 0; }	// Initialise the matrix to zero
		for(auto & pp : allPP) {
			PP[(Nseq * pp.x()) + pp.y()] = PP[(Nseq * pp.y())+pp.x()] = pp.PP(i);
		}
		// Do the divvying
		string seq;	// The column characters; needed for gaps
		for(int k = 0; k < Nseq; k++) { seq = seq + in_seq[k][i]; }
		vector <vector <int> > divvy = CCluster::Instance()->OutputClusters(PP,seq,options.threshold);
		// Sort output
		if((divvy.size() == Nseq && options.minDivvyChar > 1) && options.doFullDivvy == 0) { continue; }
		for(auto & v : divvy) {
			// Count characters for output
			int count = 0;
			for(auto x : v) { if(!IsGap(in_seq[x][i])) { count ++; } }
			// Always skip if all gaps
			if(count == 0) { continue; }
			// Skip for normal divvying or partial filtering if fewer than options.minDivvyChar characters in the column
			if(count < options.minDivvyChar && options.doFullDivvy >= 0) { continue; }
			// Do output
			colCount ++;
			for(int j = 0; j < Nseq; j++)  {
				if(IsGap(in_seq[j][i])) {	// Sort gaps
					out_seq[j] << in_seq[j][i];
				} else if(find(v.begin(),v.end(),j) != v.end()) {	// Sort  the real sequences
					out_seq[j] << in_seq[j][i];
				} else {											// Otherwise it's a divvied character
					if(options.doFullDivvy == 0) { out_seq[j] << "-"; }	// Only one column so make it look gapped
					else                         { out_seq[j] << options.divvy_char; } // Otherwise standard divvying
				}
			}
		}
	}
	// Do the output
	string outfile = fileIn + options.suffixDivvy;
	replace(outfile,".fas","");
	outfile += ".fas";
	cout << "\nDivvying complete. Outputting " << colCount << " columns to " << outfile;
	ofstream out(outfile.c_str());
	for(int i = 0; i < Nseq ; i++) {
		out << ">" << names[i] << endl;
//		out << "\nLength " << out_seq[i].str().size();
		out << out_seq[i].str() << endl;
	}
	out.close();
	if(CCluster::Instance()->Warning()) {
		cout << "\n\nWARNING: some columns had no information supporting or refuting divvying clusters";
	}
	cout << "\n\n";

}

void GetPosteriors(string File) {
	assert(Nseq > 0);
	assert(!names.empty() && !in_seq.empty());
	// If the file exists then try reading it
	if(file_exist(File)) {
		cout << "\n\tPosterior probability file <" << File << "> exists. \n\tReading file";
		vector <string> Toks;
		string tmp;
		ifstream in(File);
		// Check sequence number
		tmp = read_line(in);
		Toks = Tokenise(tmp);
		if(Toks.size() != 2) { cout << "\nNumber of sequences on first line contains multiple tokens in PP file\n"; exit(-1); }
		if(atoi(Toks[0].c_str()) != Nseq) { cout << "\nNumber of sequences in first line of PP file does not match data\n"; exit(-1); }
		if(atoi(Toks[1].c_str()) != alen) { cout << "\nLength of alignment in PP file does not match data\n"; exit(-1); }
		cout << " ... checking names" << flush;
		// Check sequence names
		for(int i = 0 ;i < Nseq; i++) {
			tmp = read_line(in);
			Toks = Tokenise(tmp);
			if(Toks.size() != 1) { cout << "\nName of sequence " << i << " has white spaces. This is not allowed\n"; exit(-1); }
			if(Toks[0] != names[i]) { cout << "\nName of sequence " << i << " does not match that in data. This is not allowed\n"; exit(-1); }
		}
		cout << " ... checking PPs" << flush;
		// Input PPs
		while(getline(in,tmp)) {
			allPP.push_back(CPostP(tmp));
		}
		in.close();
		cout << " ... success" << flush;
	}
	// Otherwise calculate the posteriors
	else {
		cout << "\n\tCalculating " << CCluster::Instance()->NoPairs() << "/" <<  ((Nseq*Nseq)-Nseq)/2 << " pairwise posteriors (min_per_split=" << options.approxNumber << "). This may take some time...\n" << flush;
		// Do the calculation
		int count = 0;
		for(auto & x : CCluster::Instance()->PairsToCalculate()) {
			ProgressSpinner(++count, CCluster::Instance()->NoPairs(),"\t");
			getSinglePosterior(x[0],x[1]);
			allPP.push_back(CPostP(x[0],x[1],zorro_posterior,alen));
		}
		// Do the output
		ofstream out(File.c_str());
		out << Nseq << " " << alen;
		for(string &n : names) { out << endl << n; }
		for(CPostP &pp : allPP) { out << endl <<pp.out(); }
		out << endl;
		out.close();
	}
}

// Wrapper function for MakeTreeNJDist to cope with weird characters in taxa names
CTree MakeTree(vector <string> n, vector <string> s) {
	vector <string> newName(n.size());
	for(int i = 0; i < n.size(); i++) { newName[i] = "seq" + int_to_string(i); }
	CTree retTree(MakeTreeNJDist(newName,s),newName);
	retTree.SetNames(n,true);
	return retTree;
}

string MakeTreeNJDist(vector <string> n, vector <string> s) {
	int no_seq = n.size();
	assert(n.size() == s.size());
	vector <double> distances(no_seq * no_seq, 0);
	for(int i = 0; i < no_seq; i++) {
		for(int j = i+1; j < no_seq; j++) {
			distances[(j*no_seq)+i] = distances[(i*no_seq)+j] = AAJCdist(GetPercentDiff(s[i],s[j]));
		}
	}
	return DoBioNJ(distances,n);
}

double AAJCdist(double p) {
	double val = 1 - ( (20.0/19.0) * p);
	if(val < DBL_EPSILON) { return BIG_NUMBER; }
	return -(19.0/20.0) * log(val);
}

double GetPercentDiff(string seq1, string seq2) {
	int diff = 0, total = 0;
	assert(seq1.size() == seq2.size());
	for(int i = 0 ; i < seq1.size(); i++) {
		if(IsGap(seq1[i]) || IsGap(seq2[i])) { continue; }
		total ++;
		if(seq1[i] != seq2[i]) { diff ++; }
	}
	if(total == 0) { return 1.0; }
	return (double) diff / (double) total;
}

void CreateZorro() {
	// Hard coded options
	JTT = 0;
	PMB = 0;
	PAM = 1;
	MATRICES = 0;
	// Other memory initialisation
	Nseq = names.size();
	alen = in_seq[0].size();
	zorro_align = (char **)(malloc(Nseq * sizeof(char *)));
	zorro_sequence = (char **)(malloc(Nseq * sizeof(char *)));
	lens = (int *)(malloc(Nseq * sizeof(int*)));
	for(int i = 0 ; i < Nseq; i++) {
		zorro_align[i] = (char*)(malloc (alen * sizeof(char*)) );
		for(int j = 0; j < in_seq[i].size(); j++) {
			zorro_align[i][j] = pep2num(in_seq[i][j]);
		}
		zorro_sequence[i] = removeGaps(zorro_align[i],alen,&lens[i]);
	}
	initHMM(alen);
	calc_prep(alen);
}

bool CheckNext(int i, int argc, char *argv[]) {
	if(i == argc) { return false; }
	if(argv[i][0] == '-') { return false; }
	return true;
}

