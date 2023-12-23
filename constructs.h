//
//    Version: $Id: constructs.h 170 2014-01-29 14:00:04Z gk $
//
/* 
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef  __CONSTRUCTS_H__
#define  __CONSTRUCTS_H__



#include <stdlib.h>
#include <sys/stat.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <pthread.h>

#include <algorithm>
#include <vector>
#include <fstream>
#include <map>


using namespace std;

enum PATTERN_MODE { 
	PATTERN_SEQ=0,
	PATTERN_SIM=1
};


const int DEBUG_NID = 91;
const int DEBUG_BID = 2270;
const int DEBUG_SID = 3024;
const int CUTOFF = 10.0;




#define TYPE_CS  1 // Conditional Stimulus 
#define TYPE_US  2
#define TYPE_NOISE  3


#define param_E_L  0.0
#define param_V_rest   param_E_L


enum  {
	RUN_TRAINING = 1,
	RUN_PRE = 2,
	RUN_TEST = 3,
};


struct LANeuron;
struct LABranch;




inline double rgen()
{
	return double(std::rand())/double(RAND_MAX);
}

/* Tracks independent synapses , synaptic tags and synaptic calcium concentration etc. */
struct LASynapse
{
	int sid;
	bool isPlastic;
	float weight, iltp, eltp, pos, calcium, stag, trigger, ltag, stdpTag ;

	float location; // 0... 1 In branch

	vector<float> weightHistory;
	vector<float> tagHistory;
	int tagTime;
	float tagMax;

	LANeuron* source_nrn;
	LANeuron* target_nrn;
	LABranch* target_branch;

	LASynapse()
	{
		Reset();
	}

	~LASynapse() 
	{ 
		Reset();
	}

	void Reset()
	{
		source_nrn = target_nrn  =0;
		target_branch = 0;
		weight = 0.9999;
		isPlastic = false;
		iltp = ltag = stag = eltp = calcium =0.0;
		pos = rgen();
		sid = -1;
		stdpTag =0;
		tagTime =-1;
	}
};



enum dend_conds {
	DEND_SUPRA = 0,
	DEND_SUB,
	DEND_LINEAR,
	DEND_MIXED,
};

/* Tracks local branch depolarization, protein production, branch strength */
struct LABranch
{
	int bid, branch_spikes ;
	LANeuron* neuron;
	float passiveWeight;
	vector<LASynapse*> synapses;
	float depol, depol2, bcalcium;
	float protein, proteinRate, strength, dreset, strengthTag, dspike, totcalc;
	int dspikeT;
	float protein1, protein2;
	vector<pair<float,float> > prpTransients;
	float turnoverPoint;
	int nlType;

	vector<float> branchStrengthHistory;	// For graphs 
	vector<float> branchProteinHistory;	// For graphs 
	vector<float> branchVoltageHistory; 		// For graphs
	vector<float> branchSpikesHistory; 		// For graphs
	vector<float> branchCalciumHistory;

	float turnoverRate;
	float dendExc, dendInh;

	LABranch() 
	{
		Reset();
	}

	~LABranch() 
	{ 
		Reset();
	}


	void Reset()
	{
		turnoverRate= 0.;
		turnoverPoint = rgen();
		bcalcium = depol =0.0;
		bid = -1;
		passiveWeight = 1.0;
		bid = -1;
		proteinRate = protein = 0.;
		strength = 1.0;
		strengthTag =0.0;
		branch_spikes = 0;
		dspike = dreset =0.0;
		depol = depol2 =0.0;
		protein1 = protein2 = 0.0;
		totcalc =0.0;
		dspikeT =0;
		nlType = DEND_SUPRA;
		dendExc= dendInh =0;
	}
};



struct LANetwork;

/* Tracks somatic potential, somatic protein availability, spiking and backpropagating */
struct LANeuron
{
	int nid, input_id, input_type;
	char type;
	float pos_x,pos_y,pos_z;
	float glx, gly;
	float V, w, crebLevel, protein, proteinRate, synScaling, totcalc;
	int crebLevelT;
	float wadapt, actvar, vreset, vspike;
	float stdp_x, branch_scaling;
	float synapticWeightsInitialSum;
	int lastSpikeT;
	float bAP;
	LANetwork* network;
	float somaExc, somaInh;



	vector<pair<float,float> > prpTransients;

	vector<float> voltageHistory; 
	vector<int>   spikeTimings; 

	vector<float>   proteinHistory; 
	vector<float>   crebHistory; 

	int total_spikes, dend_spikes;

	vector<LABranch*> branches;
	vector<LASynapse*> outgoing;
	vector<LASynapse*> incoming;

	LANeuron()
	{
		Reset();
	}

	~LANeuron() 
	{
		Reset();
	}

	void Reset()
	{
		type = ' ';
		bAP =0.0;
		w = crebLevel = protein =  wadapt = actvar = vreset= 0.0;
		nid = -1;
		V = 0.;
		crebLevelT  = -1;

		pos_x = rgen();
		pos_y = rgen();
		pos_z = rgen();
		stdp_x=0.0;
		lastSpikeT =-100;

		protein = proteinRate = 0.0;
		totcalc =0.0;

		synScaling = total_spikes =0;
		branch_scaling =0;
		input_id = -1;
		vspike =0.;
		dend_spikes =0;

		somaExc = somaInh =0;
	}

};




/* Artificial spike generators for external inputs */
struct LAInput: public LANeuron
{
	int* spikeTimes;
	int curSpike;
	int nextSpikeT;
	int totalSpikes;
	int groupIdx; // id of this neuron in the group representing a memory


	LAInput()
	{
		LANeuron();
		spikeTimes = 0;
		groupIdx =-1;
	}

	~LAInput()
	{
		Reset();
	}

	void Reset() 
	{ 
		curSpike = -1;
		delete[] spikeTimes;
		spikeTimes = 0;
		totalSpikes =0;
		nextSpikeT = -1;
	}

	

	int Program(int tstart, int duration, float freq, float randomness)
	{
		Reset();
		int total = round((float(duration) * freq)/1000.0);
		//printf("Programming %d spikes dur = %d freq=%f\n", total, duration, freq);
		if (!total) return 0;

		float period = 1000.0/freq;
		spikeTimes = new int[total];
		for (int i =0; i < total; i++)
		{
			spikeTimes[i]  =  tstart + period*i + rgen()*randomness*period;
		}

		this->curSpike =0;
		this->nextSpikeT = this->spikeTimes[this->curSpike];
		this->totalSpikes = total;
		//cout << "Total spikes " << totalSpikes << endl;
		return total;
	}


	





	int CopyShuffled(LAInput &other, float randomness)
	{
		Reset();
		this->spikeTimes = new int[other.totalSpikes];
		for (int i=0; i < other.totalSpikes; i++)
		{
			this->spikeTimes[i] = other.spikeTimes[i] += (rgen() - 0.5)*randomness;
		}

		this->totalSpikes = other.totalSpikes;
		return this->totalSpikes;
	}



};



/* This is used by wxWidgets */
struct Arr2D {
	float* data;
	int nx, ny;

	Arr2D(int nx, int ny) 
	{ 
		this->data = new float[nx*ny];
		this->nx = nx;
		this->ny = ny;
	}

	float& at(int x, int y)
	{
		return data[x*nx+y];
	}

}; 


/* used by wxWidgets */ 
class LAWindow;




/* iterator shortcuts */
typedef vector<LANeuron*> 			nrn_list;
typedef vector<vector<LANeuron*> >::iterator 	input_iter;
typedef vector<LANeuron*>::iterator 		nrn_iter;
typedef vector<LABranch*>::iterator 		branch_iter;
typedef vector<LASynapse*>::iterator 		syn_iter;
typedef vector< pair<float, float> >::iterator 	pair_iter;



/* Global structure to hold network configuration */ 
struct LANetwork
{

	public:

	vector<LASynapse*> synapses; /* List of all synapses */
	vector<LANeuron*> neurons;
	vector<LABranch*> branches;


	vector<LANeuron*> pyr_list;
	vector<LANeuron*> in_pv;
	vector<LANeuron*> in_som;
	vector<LANeuron*> da_list;

	vector<LANeuron*> noise_inputs;

	vector< vector<LANeuron*> > inputs_cs; /* List of input stimuli (each stimulus is a set of neurons) */
	vector< vector<LANeuron*> > inputs_us; /* List of input stimuli (each stimulus is a set of neurons) */

	vector<LANeuron* > inputs_binary;

	ofstream* traceFile;

	vector< vector<int> > spikesPerPattern;
	//vector< LANeuron*> us_inputs;
	vector<float> dbgNeuron;

	vector< vector<int> > spikeTimings;  // stores time of spikes during stimulation  only!
	vector< vector<int> > spikesPerStim;  // stores time of spikes during stimulation  only!
	vector< vector<float> > nrnVoltages; // ditto for voltages
	map< pair<int, int>, double> distances; // holds euclidean between neurons if needed
	static int RSEED;
	int runningMode, runningPatternNo, enableTurnover;
	int blockedLCpattern;
	float localPRPThresh, globalPRPThresh;
	float homeostasisTimeParam; // Time that it takes for synaptic scaling to be applied
	float BSPTimeParam; // Time that it takes for BSP to  be applied
	float CREBTimeParam; // Time that it takes for CREB to fall
	float spikeThreshDrop; 
	float connectivityParam, inhibitionParam, stimDurationParam; // Multipliers for doing sensitivity analysis 

	float crebDropFactor, tagsParam, proteinsParam, calciumParam, mgBlockParam;
	float blockLC;

	int setNlTypes;

	int nPyr2PV, nPyr2SOM, nPV2Pyr, nSOM2Pyr, Pyr2SOMplastic, Pyr2PVplastic;
	int nDA2Pyr, nDA2PV, nDA2SOM ;

	int weakMemId;
	int nBranchesTurnover;
	vector<int> isWeakMem;



	float injectedCurrent;
	float forceWadapt;
	char runProtocol;

	pthread_mutex_t synapses_mutex;


	ofstream mfile, vfile, sumweightsFile;

	vector< vector<int> > patterns;
	float homeostasisTime;

	int synapsesCounter;


	LAWindow* wx_window;

	FILE* spikesFile; // Save spiking info 'ere 

	Arr2D* voltageData;

	int     n_neurons, 
		n_branches_per_neuron,
		n_inputs,
		n_neurons_per_input,
		Tstimulation, // Total stimulated time
		T; // Simulation clock time 
	bool enablePlasticity; /* Is plasticity enabled? */
	bool isInterstim;
	bool disableCreb;
	bool debugMode;

	bool isRecall;

	bool localProteins, repeatedLearning, pretraining, altConnectivity, globalProteins;

	char* conditionsString;

	string datadir ;

	float branchOverlap;

	float initWeight, maxWeight, dendSpikeThresh;
	bool enablePruning, isPruningSynapses, INClustered;

	float inSomaTau , pyrSomaTau;

	int  inDendrites;

	int forceCreb=0;

	LANetwork()
	{
		blockLC=0;
		enablePruning = isPruningSynapses = false;
		synapsesCounter =0;
		T  = Tstimulation =0;
		isRecall = false;
		wx_window = NULL;
		n_neurons = n_branches_per_neuron = n_inputs = n_neurons_per_input =0;
		enablePlasticity = true;
		isInterstim = false;
		nBranchesTurnover=0;
		spikesFile = NULL;
		runningMode = RUN_TRAINING;
		runningPatternNo = 0;
		blockedLCpattern = -1;
		weakMemId = -1;
		repeatedLearning = globalProteins =  localProteins = false;
		pthread_mutex_init(&this->synapses_mutex, NULL);
		pretraining = false;
		debugMode = false;
		altConnectivity=false;
		conditionsString = NULL;
		datadir = "./";
		branchOverlap = -1.0; 
		disableCreb = false;
		homeostasisTime = 24.0;

		localPRPThresh=1.8;

		globalPRPThresh=18.0; 

		crebDropFactor = 1.0;
		tagsParam = 1.0;
		calciumParam = 1.0;
		mgBlockParam = 1.0;
		proteinsParam = 1.0;
		spikeThreshDrop = 0.0;

		homeostasisTimeParam = 1.0;
		BSPTimeParam = 1.0;
		CREBTimeParam = 1.0;
		inhibitionParam = 1.0;
		connectivityParam = 1.0;
		stimDurationParam = 1.0;
		dendSpikeThresh = 1.0;
		initWeight = 0.3;
		maxWeight = 1.0;

		traceFile = NULL;

		//turnoverHotspots =0;
		
		setNlTypes = DEND_LINEAR;
		INClustered=0;
		enableTurnover=0;

		pyrSomaTau = 30.;
		inSomaTau = 10.;

		inDendrites = 1;


		nPyr2PV = 1000;
		nPV2Pyr = 10000;

		nPyr2SOM = 1000;
		nSOM2Pyr = 2000;

		Pyr2SOMplastic = 0;
		Pyr2PVplastic = 0;

		nDA2Pyr = 0;
		nDA2PV = 2000;
		nDA2SOM = 2000;


		injectedCurrent =0;
		forceWadapt =0;
	}




	static void SetRandomSeed( int seed)
	{
		LANetwork::RSEED = seed+80; // XXX
		std::srand(seed+80);
	}



	~LANetwork() 
	{
		Cleanup();
	}

	void DoTurnover(float durationSecs );

	/* Set up the network, neurons and connections */
	void CreateFearNet(int, int , int, int);


	void Cleanup(void)
	{
		ResetSpikeCounters();

		for (nrn_iter ni = this->neurons.begin(); ni != this->neurons.end(); ni++)
			delete (*ni);

		for (branch_iter ni = this->branches.begin(); ni != this->branches.end(); ni++)
			delete (*ni);

		for (syn_iter ni = this->synapses.begin(); ni != this->synapses.end(); ni++)
			delete (*ni);

		this->neurons.clear();
		this->branches.clear();
		this->synapses.clear();

		if (this->spikesFile)
			fclose(this->spikesFile);
	}


;
	/* Connect two sets of neurons , can specify minimum  / maximum allowed distances between pairs of neurons */
	int ConnectNeurons(vector<LANeuron*> fromList, vector<LANeuron*> toList, bool isClustered, float toDistance, int nNeuronPairs, int nSynapsesPerNeuron, float weight, bool isPlastic= false, bool randomizeweight = false, float overlap =-1.0);

	void AddSynapse(LANeuron* a, LABranch* br, float weight, bool isPlastic);
	int ConnectInputs(vector<LANeuron*> fromList, vector<LANeuron*> toList, int nSynapses);

	int PurgeInputSynapses( int totalToRemove, float);
	int CreateInputSynapses( int totalToAdd);

	/* Create a set of neurons and append them to specified list */
	void CreateNeurons(int number, int n_branches_per_neuron, char type, vector<LANeuron*>* appendTo = 0, int inputId =-1, int somethingDummy = 0);

	void CalculateDistances();


	void Stimulate(int, int);

	void Stimulate2(int, int);

	void RunPatternTest();

	/* Simulate stimulus-dynamics (1msec per step) */
	void StimDynamics(int duration);

	/* Simulate inter-stimulus dynamics (protein / creb level / synapse weights changes only */
	void Interstim(int duration);

	void Begin(void);

	void CreateTags(void);

	void MemoryStats(int);

	void ResetSpikeCounters(void)
	{
		for (nrn_iter i = neurons.begin(); i != neurons.end(); i++)
		{
			(*i)->total_spikes =0;
			(*i)->dend_spikes =0;
			for (branch_iter b = (*i)->branches.begin(); b != (*i)->branches.end(); ++b)
			{
				(*b)->branch_spikes = 0;
			}
		}
	}

	void ResetCrebLevels() 
	{
		for (nrn_iter i = neurons.begin(); i != neurons.end(); i++)
			(*i)->crebLevel =0.0;
	}



	void StoreDataFiles( bool);

	void RecordInitialWeightSums()
	{
		for (nrn_iter i = neurons.begin(); i != neurons.end(); i++)
		{
			LANeuron* nrn = *i;
			nrn->synapticWeightsInitialSum =0.;
			
			for (branch_iter b = (*i)->branches.begin(); b != (*i)->branches.end(); ++b)
				for (syn_iter si = (*b)->synapses.begin(); si != (*b)->synapses.end(); ++si)
				{
					nrn->synapticWeightsInitialSum += (*si)->weight;
				}
		}
	}



	void SaveSpikeCounters(ofstream& ratesdat)
	{
		for (nrn_iter na = neurons.begin(); na != neurons.end(); ++na)
		{
			LANeuron* nrn = *na;
			ratesdat << nrn->total_spikes << " ";
		}
		ratesdat << endl;
	}

	
	void SaveSnapshot(char*);

	bool HasCondition(char* cond)
	{
		if (this->conditionsString && strstr(this->conditionsString, cond))
			return true;
		return false;
	}

	bool SaveSynapseState(char* filename)
	{
		ofstream synstatedat(filename);

		for (syn_iter si =this->synapses.begin(); si != this->synapses.end(); ++si)
		{
			LASynapse* s = *si;
			if (s->isPlastic)
			{
			synstatedat << s->sid<<" " 
					<< s->target_branch->bid<<" "
					<< s->target_nrn->nid  << " "
					<< s->source_nrn->nid <<" " 
					<< s->source_nrn->input_id<< " "
					<< s->target_branch->strength  << " "
					<< s->weight << " " 
					<<endl;


			}
		}
		return true;
	}

	void ReportSumWeights()
	{

		float consolidated[this->n_inputs] ;

		for (int i=0; i < this->n_inputs; i++) consolidated[i] =0;

		if (!this->sumweightsFile.is_open())
		{
			this->sumweightsFile.open((this->datadir + "/sum-weights.txt").c_str(), std::ofstream::out );
		}

		int totalPot =0;
		for (syn_iter si = this->synapses.begin(); si != this->synapses.end(); ++si)
		{
			LASynapse* s =*si;
			if ( s->source_nrn->input_id >=0 && s->target_nrn->type == 'P')
			{
				if (s->weight > 0.7)
					totalPot ++;
				consolidated[s->source_nrn->input_id] += s->weight;
			}
		}
		
		for (int i=0; i < this->n_inputs; i++)
		{
			//cout << " ["<< i <<"]-> "<< consolidated[i]<< endl;

			this->sumweightsFile << consolidated[i] <<  " " ;
			cout << consolidated[i] <<  " " ;
		}
		this->sumweightsFile << endl;
		cout << endl;

		cout <<"Total psyn: "<<totalPot<<endl;
	}

	void SetDataDir(string dir)
	{
		this->datadir = dir;
		mkdir(this->datadir.c_str(), 0755);

	}

	void PrintSynapsesSnapshot(string outfile)
	{
		ofstream  fout(outfile.c_str());

		for (nrn_iter ni = this->pyr_list.begin(); ni != this->pyr_list.end(); ni++)
		{
			LANeuron*  nrn = *ni;
			for (branch_iter bi = nrn->branches.begin(); bi != nrn->branches.end(); ++bi)
			{
				LABranch* b = *bi;
				for (syn_iter si = b->synapses.begin(); si != b->synapses.end(); ++si)
				{
					LASynapse* s =*si;
					LAInput* src = (LAInput*)s->source_nrn;
					if (src->input_id >=0) 
					{
						fout  << src->input_id << " " << src->groupIdx << " " << s->target_branch->bid << " " << s->target_nrn->nid << " " << s->weight << " "  << endl;
					}
					
				}
			}
		}
	}


	void SaveCalcs()
	{
		ofstream ooo("./data/calc.dat");
		for (syn_iter si = this->synapses.begin(); si!= this->synapses.end(); si++)
		{
			LASynapse* s = *si;
			if (s->source_nrn->input_id ==0)
			{
			//	ooo << s->calcium << endl;
			}

		}
		ooo.close();
	}


	void RunTests();
	
};





#endif
