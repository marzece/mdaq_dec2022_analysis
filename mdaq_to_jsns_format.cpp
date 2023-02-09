// Compile: 'g++ -o mdaq_to_jsns_format mdaq_to_jsns_format.cpp `root-config --cflags --libs`
// Run: './mdaq_to_jsns_format INPUT_FILENAME OUTPUT_FILENAME'
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <assert.h>
#include <arpa/inet.h>

// This is the amount of space that will be allocated to hold the
// data for a single event. Corresponds to 1MB per channel
#define DATA_MAX_SIZE (1024*1024*128)

#define htonll(x) ((((uint64_t)htonl(x)) << 32) + htonl((x) >> 32))
#define ntohll(x) htonll(x)


// Below #defines set the number of channels this script will look for in each
// event. If the data does not have all 128 channels of data, these should be modified.
#define NUM_CHANNELS_PER_MOD  16
#define MAX_NUM_MOD 32
#define FONTUS_DEVICE_MASK 0x1


#define EB_HEADER_SIZE 16
typedef struct EB_Header {
    uint32_t trig_number;
    uint32_t status;
    uint64_t device_mask;
} EB_Header;

#define XEM_HEADER_SIZE  20
typedef struct XEM_Header {
    uint32_t magic_number;
    uint32_t trig_number;
    uint64_t clock;
    uint16_t length;
    uint8_t device_id;
    uint8_t crc;
} XEM_Header;

#define FONTUS_HEADER_SIZE  52
typedef struct FONTUS_Header {
    uint32_t magic_number;
    uint32_t trig_number;
    uint64_t timestamp;
    uint16_t trig_length; 
    uint8_t device_id; 
    uint8_t trig_flag; 
    uint32_t self_trigger_word; 
    uint64_t beam_time; 
    uint64_t led_trigger_time; 
    uint64_t ct_time; 
    uint32_t crc; 
} FONTUS_Header;

typedef struct ROOTOutData {
    // FONTUS info
    uint32_t device_mask;
    uint32_t trigger_word;
    uint32_t self_trigger_word;
    double fontus_timetag;
    double led_timetag;
    double kicker_timetag;
    double ct_timetag;
    uint32_t nmod;
    // CERES info
    uint32_t nsample;
    uint32_t nchan;
    double* timetag;
    uint16_t* fadc;
} ROOTOutData;

int main(int argc, char** argv) {
    if(argc != 3) {
        printf("Must provide an input and output filename\n");
        return 0;
    }
    char* filename = argv[1];
    FILE* fin = fopen(filename, "rb");
    if(!fin) {
        printf("Could not open file \"%s\":\n", filename);
        perror(NULL);
        return -1;
    }
    TFile* fout = new TFile(argv[2], "RECREATE");

    ROOTOutData out_data;
    ROOTOutData temp_out_data;
    memset(&temp_out_data, 0, sizeof(ROOTOutData));
    memset(&out_data, 0, sizeof(ROOTOutData));

    out_data.fadc = (uint16_t*)malloc(DATA_MAX_SIZE);
    out_data.timetag = (double*)malloc(sizeof(double)*MAX_NUM_MOD);

    // These two structs should share their dynamic buffers b/c just trust me
    temp_out_data.fadc = out_data.fadc;
    temp_out_data.timetag = out_data.timetag;

    TTree * tree = new TTree("tree", "tree");
    TBranch* trigger_word_branch = tree->Branch("trigger_word", &out_data.trigger_word, "trigger_word/i");
    TBranch* self_trigger_word_branch = tree->Branch("self_trigger_word", &out_data.self_trigger_word, "self_trigger_word/i");
    TBranch* device_mask_branch = tree->Branch("device_mask", &out_data.device_mask, "device_mask/i");
    TBranch* led_timestamp_branch = tree->Branch("led_timestamp", &out_data.led_timetag, "led_timestamp/D");
    TBranch* kicker_timestamp_branch = tree->Branch("kicker_timestamp", &out_data.kicker_timetag, "kicker_timestamp/D");
    TBranch* ct_timestamp_branch = tree->Branch("ct_timestamp", &out_data.ct_timetag, "ct_timestamp/D");
    TBranch* nmod_branch = tree->Branch("nmod", &out_data.nmod, "nmod/i");
    TBranch* nsample_branch = tree->Branch("nsample", &out_data.nsample, "nsample/i");
    TBranch* timetag_branch = tree->Branch("TimeTag", out_data.timetag, "TimeTag[nmod]/D");
    TBranch* nchan_branch = tree->Branch("nchan", &out_data.nchan, "nchan/i");
    TBranch* fadc_branch = tree->Branch("FADC", out_data.fadc, "FADC[nsample][128]/s");

    int i, j, k;
    unsigned char buffer[1024];
    FONTUS_Header trig_header;
    memset(&trig_header, 0, sizeof(trig_header));
    XEM_Header headers[16];
    long data_start_locations[16];
    EB_Header eb_header;


    // Get the data file size
    fseek(fin, 0, SEEK_END);
    long file_size = ftell(fin);
    rewind(fin);

    // Read events
    // Start with the event builder header
    int loop=0;
    while(1) {
        printf("%i\n", loop++);
        if(fread(buffer, EB_HEADER_SIZE, 1, fin) !=1) {
            //Error
            break;
        }
        eb_header.device_mask = htonll(*((uint64_t*)(buffer + 8)));
        // We'll only look at the bottom 32-bits
        temp_out_data.device_mask = (uint32_t)(eb_header.device_mask & 0xFFFFFFFF);
        if(temp_out_data.device_mask != eb_header.device_mask) {
            printf("Non 32-bit device_mask encoutered. I can't handle this!\n");
            break;
        }


        temp_out_data.nmod = 0;
        // Count the number of bits set in the device mask to get the number of
        // modules
        for(uint32_t temp=(temp_out_data.device_mask & ~FONTUS_DEVICE_MASK); temp; temp=temp>>1) {
            temp_out_data.nmod += temp & 0x1;
        }

        temp_out_data.nmod = 8;// TODO
        temp_out_data.nchan = temp_out_data.nmod*NUM_CHANNELS_PER_MOD;


        // Then the FONTUS header
        if(temp_out_data.device_mask & FONTUS_DEVICE_MASK) {
            if(fread(buffer, FONTUS_HEADER_SIZE, 1, fin) !=1) {
                //Error
                break;
            }
            trig_header.magic_number = ntohl(*((uint32_t*)buffer));
            trig_header.trig_number = ntohl(*((uint32_t*) (buffer+4)));
            trig_header.timestamp = ntohll(*((uint64_t*) (buffer+8)));
            trig_header.trig_length = ntohs(*((uint16_t*) (buffer+16)));
            trig_header.device_id = *((uint8_t*) (buffer+18)); 
            trig_header.trig_flag = *((uint8_t*) (buffer+19)); 
            trig_header.self_trigger_word = ntohl(*((uint32_t*) (buffer+20)));
            trig_header.beam_time = *((uint64_t*) (buffer+24)); 
            trig_header.led_trigger_time = *((uint64_t*) (buffer+32)); 
            trig_header.ct_time = *((uint64_t*) (buffer+40)); 
            trig_header.crc = *((uint32_t*) (buffer+48));

            temp_out_data.fontus_timetag = trig_header.timestamp;
            temp_out_data.trigger_word = trig_header.trig_flag;
            temp_out_data.self_trigger_word = trig_header.self_trigger_word;
            temp_out_data.led_timetag = trig_header.led_trigger_time;
            temp_out_data.kicker_timetag = trig_header.beam_time;
            temp_out_data.ct_timetag = trig_header.ct_time;

            // See if we have to do the look back thing
            if(loop!=0) {
                double earliest_time = 0;
                double latest_time = 0;
                for(int i=0; i<temp_out_data.nmod; i++) {
                    double this_time = temp_out_data.timetag[i];

                    if(i==0 || this_time < earliest_time) {
                        earliest_time = this_time;
                    }
                    if(i==0 || this_time > latest_time) {
                        latest_time = this_time;
                    }
                }
                // Each sample represents 2ns, but timestamp the clock runs at
                // 250MHz, which corresponds to 4ns per tick.
                latest_time += temp_out_data.nsample/2.;

                if(temp_out_data.led_timetag < latest_time && temp_out_data.led_timetag > earliest_time) {
                    temp_out_data.led_timetag = temp_out_data.led_timetag;
                    temp_out_data.led_timetag = 0;
                }
                if(temp_out_data.ct_timetag < latest_time && temp_out_data.ct_timetag > earliest_time) {
                    temp_out_data.ct_timetag = temp_out_data.ct_timetag;
                    temp_out_data.ct_timetag = 0;
                }
                if(temp_out_data.kicker_timetag < latest_time && temp_out_data.kicker_timetag > earliest_time) {
                    temp_out_data.kicker_timetag = temp_out_data.kicker_timetag;
                    temp_out_data.kicker_timetag = 0;
                }
            }
        }

        if(loop!=0) {
            tree->Fill();
        }

        out_data = temp_out_data;
        // Now go to each CERES XEM header and collect the headers
        for(i=0; i<out_data.nmod; i++) {
            if(fread(buffer, XEM_HEADER_SIZE, 1, fin) != 1) {
                //Error
                goto DONE;
            }
            data_start_locations[i] = ftell(fin);

            headers[i].magic_number = ntohl(*((uint32_t*)buffer));
            headers[i].trig_number = ntohl(*((uint32_t*) (buffer+4)));
            headers[i].clock = ntohll(*((uint64_t*) (buffer+8)));
            headers[i].length = ntohs(*((uint16_t*) (buffer+16)));
            headers[i].device_id = *((uint8_t*) (buffer+18)); 
            headers[i].crc = *((uint8_t*) (buffer+19));

            assert(headers[i].magic_number == 0xFFFFFFFF);

            out_data.timetag[i] = (double)headers[i].clock;
            // If this is not the last module, step forward to the next XEM header
            if(i != out_data.nmod-1) {
                if(fseek(fin, (headers[i].length+2)*4*NUM_CHANNELS_PER_MOD , SEEK_CUR)) {
                    // Error
                    goto DONE;
                }
            }
        }
        // Now need to find the smallest waveform length for this event;
        uint32_t nsample = headers[0].length;
        for(i=0; i<out_data.nmod; i++) {
            if(nsample > headers[i].length) {
                nsample = headers[i].length;
            }
        }

        // The data sample count is for number of 32-bit words.
        // But each of those 32-bits is composed of 2 16-bit samples.
        out_data.nsample = nsample*2;

        for(i=0; i<out_data.nmod; i++) {
            // Go to the start of the actual ADC data
            if(fseek(fin, data_start_locations[i], SEEK_SET)) {
                //Error
                goto DONE;
            }
            for(j=0; j<NUM_CHANNELS_PER_MOD; j++) {
                // Go to the start of data for the next waveform (if there is one)
                    if(fseek(fin, data_start_locations[i] + j*(headers[i].length +2)*sizeof(uint32_t), SEEK_SET)) {
                        // Error
                        goto DONE;
                    }
                // First skip over the channel header
                if(fseek(fin, sizeof(uint32_t), SEEK_CUR) != 0) {
                    // Error
                    goto DONE;
                }

                uint16_t* data_loc = out_data.fadc + (i*NUM_CHANNELS_PER_MOD + j)*nsample;
                if(fread(data_loc, sizeof(uint16_t), nsample, fin) != nsample) {
                    // Error
                    goto DONE;
                }
                for(k=0; k<nsample;k++) {
                    data_loc[k] = ntohs(data_loc[k]);
                }
            }
        }
        //tree->Fill();

        long next_event = data_start_locations[out_data.nmod-1] + (headers[out_data.nmod-1].length+2)*sizeof(uint32_t)*NUM_CHANNELS_PER_MOD;
        if(fseek(fin, next_event, SEEK_SET)) {
            goto DONE;
        }
    }
    // Make sure to get the final event
    tree->Fill();


DONE:

    tree->Write();
    fout->Close();

    free(out_data.fadc);
    free(out_data.timetag);

    return 0;
}
