#!/bin/bash

DATA_FOLDER="/n/groups/marks/projects/indels_human/data/indels_raw_DMS"

# # Russ 2020 
TARGET_SEQ="MTSENPLLALREKISALDEKLLALLAERRELAVEVGKAKLLSHRPVRDIDRERDLLERLITLGKAHHLDAHYITRLFQLIIEDSVLTQQALLQQHLNKINPHSARIAFLGPKGSYSHLAARQYAARHFEQFIESGCAKFADIFNQVETGQADYAVVPIENTSSGAINDVYDLLQHTSLSIVGEMTLTIDHCLLVSGTTDLSTINTVYSHPQPFQQCSKFLNRYPHWKIEYTESTSAAMEKVAQAKSPHVAALGSEAGGTLYGLQVLERIEANQRQNFTRFVVLARKAINVSDQVPAKTTLLMATGQQAGALVEALLVLRNHSLIMTRLESRPIHGNPWEEMFYLDIQANLESAEMQKALKELGEITRSMKVLGCYPSENVVPVDPT"
SEQ_SUFFIX="LNKINPHSARIAFLGPKGSYSHLAARQYAARHFEQFIESGCAKFADIFNQVETGQADYAVVPIENTSSGAINDVYDLLQHTSLSIVGEMTLTIDHCLLVSGTTDLSTINTVYSHPQPFQQCSKFLNRYPHWKIEYTESTSAAMEKVAQAKSPHVAALGSEAGGTLYGLQVLERIEANQRQNFTRFVVLARKAINVSDQVPAKTTLLMATGQQAGALVEALLVLRNHSLIMTRLESRPIHGNPWEEMFYLDIQANLESAEMQKALKELGEITRSMKVLGCYPSENVVPVDPT"
INPUT_FILE="$DATA_FOLDER/B1LPA6_ECOSM_Russ_2020.csv" 
OUTPUT_FILE="../../processed_data/DMS/B1LPA6_ECOSM_Russ_2020.csv"
MUTANT_COLUMN="mutant" 
PHENOTYPE_NAME="activity" 
DIRECTIONALITY=1

python step1_clean_mutational_scan.py --input_file $INPUT_FILE --output_file $OUTPUT_FILE --DMS_mutant_column $MUTANT_COLUMN \
--DMS_phenotype_name $PHENOTYPE_NAME --DMS_directionality $DIRECTIONALITY --target_seq $TARGET_SEQ --seq_suffix $SEQ_SUFFIX

# Gonzalez 2019 
TARGET_SEQ="MSIQHFRVALIPFFAAFCLPVFAHPETLVKVKDAEDQLGARVGYIELDLNSGKILESFRPEERFPMMSTFKVLLCGAVLSRVDAGQEQLGRRIHYSQNDLVEYSPVTEKHLTDGMTVRELCSAAITMSDNTAANLLLTTIGGPKELTAFLHNMGDHVTRLDRWEPELNEAIPNDERDTTMPAAMATTLRKLLTGELLTLASRQQLIDWMEADKVAGPLLRSALPAGWFIADKSGAGERGSRGIIAALGPDGKPSRIVVIYTTGSQATMDERNRQIAEIGASLIKHW"
INPUT_FILE="$DATA_FOLDER/BLAT_ECOLX_Gonzalez_indels_2019.csv" 
OUTPUT_FILE="../../processed_data/DMS/BLAT_ECOLX_Gonzalez_indels_2019.csv"
MUTANT_COLUMN="sequence" 
PHENOTYPE_NAME="DMS_score" 
DIRECTIONALITY=1 

python step1_clean_mutational_scan.py --input_file $INPUT_FILE --output_file $OUTPUT_FILE --DMS_mutant_column $MUTANT_COLUMN \
--DMS_phenotype_name $PHENOTYPE_NAME --DMS_directionality $DIRECTIONALITY --target_seq $TARGET_SEQ

# Sinai 2021
TARGET_SEQ="MAADGYLPDWLEDTLSEGIRQWWKLKPGPPPPKPAERHKDDSRGLVLPGYKYLGPFNGLDKGEPVNEADAAALEHDKAYDRQLDSGDNPYLKYNHADAEFQERLKEDTSFGGNLGRAVFQAKKRVLEPLGLVEEPVKTAPGKKRPVEHSPVEPDSSSGTGKAGQQPARKRLNFGQTGDADSVPDPQPLGQPPAAPSGLGTNTMATGSGAPMADNNEGADGVGNSSGNWHCDSTWMGDRVITTSTRTWALPTYNNHLYKQISSQSGASNDNHYFGYSTPWGYFDFNRFHCHFSPRDWQRLINNNWGFRPKRLNFKLFNIQVKEVTQNDGTTTIANNLTSTVQVFTDSEYQLPYVLGSAHQGCLPPFPADVFMVPQYGYLTLNNGSQAVGRSSFYCLEYFPSQMLRTGNNFTFSYTFEDVPFHSSYAHSQSLDRLMNPLIDQYLYYLSRTNTPSGTTTQSRLQFSQAGASDIRDQSRNWLPGPCYRQQRVSKTSADNNNSEYSWTGATKYHLNGRDSLVNPGPAMASHKDDEEKFFPQSGVLIFGKQGSEKTNVDIEKVMITDEEEIRTTNPVATEQYGSVSTNLQRGNRQAATADVNTQGVLPGMVWQDRDVYLQGPIWAKIPHTDGHFHPSPLMGGFGLKHPPPQILIKNTPVPANPSTTFSAAKFASFITQYSTGQVSVEIEWELQKENSKRWNPEIQYTSNYNKSVNVDFTVDTNGVYSEPRPIGTRYLTRNL" 
INPUT_FILE="$DATA_FOLDER/CAPSD_AAV2S_Sinai_indels_2021.csv" 
OUTPUT_FILE="../../processed_data/DMS/CAPSD_AAV2S_Sinai_indels_2021.csv"
MUTANT_COLUMN="full_sequence"
PHENOTYPE_NAME="viral_selection" 
DIRECTIONALITY=1

python step1_clean_mutational_scan.py --input_file $INPUT_FILE --output_file $OUTPUT_FILE --DMS_mutant_column $MUTANT_COLUMN \
--DMS_phenotype_name $PHENOTYPE_NAME --DMS_directionality $DIRECTIONALITY --target_seq $TARGET_SEQ

# Pokusaeva 2019 
TARGET_SEQ="MTEQKALVKRITNETKIQIAISLKGGPLAIEHSIFPEKEAEAVAEQATQSQVINVHTGIGFLDHMIHALAKHSGWSLIVECIGDLHIDDHHTTEDCGIALGQAFKEALGAVRGVKRFGSGFAPLDEALSRAVVDLSNRPYAVVELGLQREKVGDLSCEMIPHFLESFAEASRITLHVDCLRGKNDHHRSESAFKALAVAIREATSPNGTNDVPSTKGVLM"
INPUT_FILE="$DATA_FOLDER/HIS7_YEAST_Pokusaeva_indels_2019.csv" 
OUTPUT_FILE="../../processed_data/DMS/HIS7_YEAST_Pokusaeva_indels_2019.csv"
MUTANT_COLUMN="sequence"
PHENOTYPE_NAME="DMS_score"
DIRECTIONALITY=1 

python step1_clean_mutational_scan.py --input_file $INPUT_FILE --output_file $OUTPUT_FILE --DMS_mutant_column $MUTANT_COLUMN \
--DMS_phenotype_name $PHENOTYPE_NAME --target_seq $TARGET_SEQ

# Mighell 2018 
TARGET_SEQ="MTAIIKEIVSRNKRRYQEDGFDLDLTYIYPNIIAMGFPAERLEGVYRNNIDDVVRFLDSKHKNHYKIYNLCAERHYDTAKFNCRVAQYPFEDHNPPQLELIKPFCEDLDQWLSEDDNHVAAIHCKAGKGRTGVMICAYLLHRGKFLKAQEALDFYGEVRTRDKKGVTIPSQRRYVYYYSYLLKNHLDYRPVALLFHKMMFETIPMFSGGTCNPQFVVCQLKVKIYSSNSGPTRREDKFMYFEFPQPLPVCGDIKVEFFHKQNKMLKKDKMFHFWVNTFFIPGPEETSEKVENGSLCDQEIDSICSIERADNDKEYLVLTLTKNDLDKANKDKANRYFSPNFKVKLYFTKTVEEPSNPEASSSTSVTPDVSDNEPDHYRYSDTTDSDPENEPFDEDQHTQITKV"
INPUT_FILE="$DATA_FOLDER/PTEN_HUMAN_Mighell_deletions_2018.csv"
OUTPUT_FILE="../../processed_data/DMS/PTEN_HUMAN_Mighell_deletions_2018.csv"
MUTANT_COLUMN="sequence"
PHENOTYPE_NAME="DMS_score"
DIRECTIONALITY=1 

python step1_clean_mutational_scan.py --input_file $INPUT_FILE --output_file $OUTPUT_FILE --DMS_mutant_column $MUTANT_COLUMN \
--DMS_phenotype_name $PHENOTYPE_NAME --target_seq $TARGET_SEQ

# Kotler 2018 
TARGET_SEQ="MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD"
INPUT_FILE="$DATA_FOLDER/P53_HUMAN_Kotler_deletions_2018.csv" 
OUTPUT_FILE="../../processed_data/DMS/P53_HUMAN_Kotler_deletions_2018.csv"
MUTANT_COLUMN="sequence"
PHENOTYPE_NAME="RFS_H1299" 
DIRECTIONALITY=-1

python step1_clean_mutational_scan.py --input_file $INPUT_FILE --output_file $OUTPUT_FILE --DMS_mutant_column $MUTANT_COLUMN \
--DMS_phenotype_name $PHENOTYPE_NAME --DMS_directionality $DIRECTIONALITY --target_seq $TARGET_SEQ

TARGET_SEQ="MGSVRTNRYSIVSSEEDGMKLATMAVANGFGNGKSKVHTRQQCRSRFVKKDGHCNVQFINVGEKGQRYLADIFTTCVDIRWRWMLVIFCLAFVLSWLFFGCVFWLIALLHGDLDTSKVSKACVSEVNSFTAAFLFSIETQTTIGYGFRCVTDECPIAVFMVVFQSIVGCIIDAFIIGAVMAKMAKPKKRNETLVFSHNAVIAMRDGKLCLMWRVGNLRKSHLVEAHVRAQLLKSRITSEGEYIPLDQIDINVGFDSGIDRIFLVSPITIVHEIDEDSPLYDLSKQDIDNADFEIVVILEGMVEATAMTTQCRSSYLANEILWGHRYEPVLFEEKHYYKVDYSRFHKTYEVPNTPLCSARDLAEKKYILSNANSFCYENEVALTSKEEEEDSENGVPESTSTDSPPGIDLHNQASVPLEPRPLRRESEI"
INPUT_FILE="$DATA_FOLDER/KCNJ2_MOUSE_Macdonald_2022_indels_noflag.csv" 
OUTPUT_FILE="../../processed_data/DMS/KCNJ2_MOUSE_Macdonald_2022_indels_noflag.csv"
MUTANT_COLUMN="mutated_sequence_no_flag"
PHENOTYPE_NAME="score" 
DIRECTIONALITY=1

python step1_clean_mutational_scan.py --input_file $INPUT_FILE --output_file $OUTPUT_FILE --DMS_mutant_column $MUTANT_COLUMN \
--DMS_phenotype_name $PHENOTYPE_NAME --DMS_directionality $DIRECTIONALITY --target_seq $TARGET_SEQ

TARGET_SEQ="MKNCLKMKNLLPALTITMAMSAVMALVVTPNAYASKWDEKMTPEQVEATLDKKFAEGNYSPKGADSCLMCHKKSEKVMDLFKGVHGAIDSSKSPMAGLQCEACHGPLGQHNKGGNEPMITFGKQSTLSADKQNSVCMSCHQDDKRMSWNGGHHDNADVACASCHQVHVAKDPVLSKNTEMEVCTSCHTKQKADMNKRSSHPLKWAQMTCSDCHNPHGSMTDSDLNKPSVNDTCYSCHAEKRGPKLWEHAPVTENCVTCHNPHGSVNDGMLKTRAPQLCQQCHASDGHASNAYLGNTGLGSNVGDNAFTGGRSCLNCHSQVHGSNHPSGKLLQR"
INPUT_FILE="$DATA_FOLDER/Q8EG35_SHEON_Campbell_2022.csv" 
OUTPUT_FILE="../../processed_data/DMS/Q8EG35_SHEON_Campbell_2022.csv"
MUTANT_COLUMN="mutated_sequence"
PHENOTYPE_NAME="selected_avg" 
DIRECTIONALITY=1

python step1_clean_mutational_scan.py --input_file $INPUT_FILE --output_file $OUTPUT_FILE --DMS_mutant_column $MUTANT_COLUMN \
--DMS_phenotype_name $PHENOTYPE_NAME --DMS_directionality $DIRECTIONALITY --target_seq $TARGET_SEQ

OUTPUT_DIR="../../processed_data/DMS"

python step2_make_mapfile.py -input_folder=$OUTPUT_DIR -mapping_file=$OUTPUT_DIR/mapfiles/DMS_mapping.csv --alignment_folder="../../alignments/focus_column_only/DMS"

# chunking CAPSD DMS 

INPUT_FILE="../../processed_data/DMS/CAPSD_AAV2S_Sinai_indels_2021.csv"
OUTPUT_DIR="../../processed_data/DMS_chunked/CAPSD_AAV2S_Sinai_indels_2021"
NUM_CHUNKS=100
DMS_ID="CAPSD_AAV2S_Sinai_indels_2021"
ORIGINAL_MAPPING_FILE="../../processed_data/DMS/mapfiles/DMS_mapping.csv"

python chunk_mutational_scan.py --input_file $INPUT_FILE --output_dir $OUTPUT_DIR --num_chunks $NUM_CHUNKS \
--target_sequence_id $DMS_ID --original_mapping_file $ORIGINAL_MAPPING_FILE