//
// Version: $Id: lamodel.cpp 172 2014-02-12 10:06:07Z gk $
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


#include "constructs.h"
#include <iostream>
#include <cstring>
#include <string>
#include <unistd.h>
#include <getopt.h>
#include <iterator>




// Preallocates spikes in a list  of neurons
inline void program_input(nrn_list lst, int tstart, int duration, float freq, float randomness, int limitActive = -1)
{
	int skip = 0;
	if (limitActive == -999) skip = 3;
	for (nrn_iter n = lst.begin(); n != lst.end(); ++n)
	{
		if (skip>0)
		{
			skip--;
			continue;
		}

		LAInput* in = (LAInput*)(*n);
		in->Program(tstart, duration, freq, randomness);
	}
}





void RunPattern(LANetwork* net, vector<int>& pattern, float hifreq,  float lowfreq, int duration, int limitActive )
{

	int pad = (4000 - duration)/2;
	//printf("Running pattern: [ ");


	for (size_t j=0; j < pattern.size(); j++)
	{
		if (pattern[j]>0)
			program_input(net->inputs_cs[j], pad, duration, hifreq, 0.5, limitActive);
		else
			program_input(net->inputs_cs[j], pad, duration, lowfreq, 0.5, limitActive);

	}


	/*
	 * blockedLCpattern == -1: all patterns blocked
	 * blockedLCpattern == -2: NO patterns blocked
	 * blockedLCpattern == 0:  LC is blocked for first pattern
	 * blockedLCpattern == 1:  LC is blocked for second pattern
	 */
	if (net->blockedLCpattern == -1 || net->runningPatternNo  == net->blockedLCpattern)
	{

		net->spikeThreshDrop = 0.;


		// Default value is 1.0 . value<1.0 increases firing rate (lower adaptation), value >1.0 decreases firing rate (higher adaptation)
		net->crebDropFactor = 1.0;

		// Default is 1.0. Multiplies the amplitude of synaptic tags (both LTP and LTD tags)
		net->tagsParam = 1.0; 

		// Default is 1.0. Multiplies the amount of calcium that is entered after a presynaptic spike is combined  with postsynaptic depolarization
		net->calciumParam = 1.0;

		// Default is 1. Multiples the level of plasticity-related proteins during Interstimulus (consolidation)
		net->proteinsParam = 1.0; 

		// Default is 1. Multiplies the threshold for NMDA activation (Mg block). <1: lower threshold , >1 : Higher threshold
		net->mgBlockParam = 1.0; 

		// Set activity level of DA neurons to low frequency (0.1)
	}
	else
	{
		net->crebDropFactor = 1.0;
		net->tagsParam = 1.0;
		net->calciumParam = 1.0;
		net->proteinsParam = 1.0;
		net->mgBlockParam = 1.0;

		// High activity level of DA neurons // 30Hz
	}

	program_input(net->da_list, pad, duration, 80., 0.5, -1);


	net->ResetSpikeCounters();
	net->StimDynamics(duration+pad+pad);
	
	

	int tActive =0;
	int tSpikes =0;

	for (nrn_iter ni = net->pyr_list.begin(); ni != net->pyr_list.end(); ni++)
	{
		LANeuron* n = *ni;
		if (float(n->total_spikes) /4000.> CUTOFF )
			tActive++;
		tSpikes += n->total_spikes;
		cout << n->total_spikes <<"," ;
	}
	cout << endl;
	cout << "INpv=";
	for (nrn_iter ni = net->in_pv.begin(); ni != net->in_pv.end(); ni++)
	{
		LANeuron* n = *ni;
		cout << n->total_spikes <<"," ;
	}
	cout << endl;
	/*
	
	printf("Active: %d (%.2f%%), avgF: %.2f\n", tActive, 100.0*float(tActive)/float(net->pyr_list.size()), tSpikes/(float(net->pyr_list.size())*2));
	*/
}


void Stimulate(LANetwork* net, int mem, int hifreq)
{
	int pad=100;
	int duration=3800;

	program_input(net->inputs_cs[mem], 	      pad, duration, hifreq, 0.5 );
	program_input(net->inputs_cs[mem==0 ? 1 : 0], pad, duration, 0,   0.5 );

	cout<< "Running mem #" << mem << " ff "<<hifreq<< endl;
	net->ResetSpikeCounters();
	net->StimDynamics(duration+pad+pad);

	for (nrn_iter ni = net->pyr_list.begin(); ni != net->pyr_list.end(); ni++)
	{
		LANeuron* nrn = *ni;
		cout << nrn->total_spikes <<",";
	}
	cout<<endl;
	cout << "INpv=";
	for (nrn_iter ni = net->in_pv.begin(); ni != net->in_pv.end(); ni++)
	{
		LANeuron* n = *ni;
		cout << n->total_spikes <<"," ;
	}
	cout << endl;
	cout << "INsom=";
	for (nrn_iter ni = net->in_som.begin(); ni != net->in_som.end(); ni++)
	{
		LANeuron* n = *ni;
		cout << n->total_spikes <<"," ;
	}
	cout << endl;


}


void ResetCrebLevels(LANetwork* net)
{

	for (nrn_iter ni = net->pyr_list.begin(); ni != net->pyr_list.end(); ni++)
	{
		LANeuron* nrn = *ni;
		nrn->crebLevel =0.;
	}
}



// The main protocol for storing / homeostasis/recalling memories sequentially
void RunStoreTest(LANetwork* net, int n_patterns, int activePerPattern, int interstimMins)
{
	printf("Running net with %d pyr. neurons, %d branches, %d synapses [%s,%s] \n", (int)net->pyr_list.size(),  (int)net->branches.size(),  (int)net->synapses.size(), net->localProteins ? "Local" : "Global", net->disableCreb ? "-CREB" : "+CREB" );


	cout << "Running Training .. " << endl;
	net->enablePlasticity = true; // Enable plasticity during training


	net->spikesPerPattern.resize(n_patterns);
	for (int i=0; i < 2; i++)
		net->spikesPerPattern[i].resize(net->neurons.size());


	float stimFreq = 80.;
	if (net->nDA2Pyr) stimFreq = 75.;

	Stimulate(net, 0, stimFreq);
	for (nrn_iter na = net->neurons.begin(); na != net->neurons.end(); ++na)
	{
		LANeuron* nrn = *na;
		net->spikesPerPattern[0][nrn->nid] = nrn->total_spikes;
	}


	cout<< "Delay "<< interstimMins << "minutes" << endl;
	net->Interstim(interstimMins*60);

	

	int crebTotals = 0;
	if (!net->disableCreb)
	{
		int toReplace = 0;
		for (nrn_iter ni = net->pyr_list.begin(); ni != net->pyr_list.end(); ni++)
		{
			LANeuron* nrn = *ni;
			if (toReplace >0 && nrn->crebLevel  >0)
			{
				for (nrn_iter ni2 = net->pyr_list.begin(); ni2 != net->pyr_list.end(); ni2++)
				{
					LANeuron* nrn2 = *ni2;
					if (nrn2->crebLevel<=0)
					{
						nrn->crebLevel =0;
						nrn2->crebLevel =0;
						toReplace--;
						break;
					}
				}
			}

			if (nrn->crebLevel >0) crebTotals++;
		}
	}

	printf("CREB UP IN %d / %d\n\n", crebTotals, (int)net->pyr_list.size());
	
	Stimulate(net, 1, stimFreq);

	for (nrn_iter na = net->neurons.begin(); na != net->neurons.end(); ++na)
	{
		LANeuron* nrn = *na;
		net->spikesPerPattern[1][nrn->nid] = nrn->total_spikes;
	}




	/*

	cout << "Delay 24hours..." << endl;
	//ResetCrebLevels(net);
	//net->Interstim((int)(24*3600.));

	net->enablePlasticity = false; // No plasticity during recall

	// Recall the memories now
	cout << "Recalling memories " << endl;

	net->enablePlasticity = false;
	ResetCrebLevels(net);

	net->spikesPerPattern.resize(n_patterns);
	for (int i=0; i < 2; i++)
		net->spikesPerPattern[i].resize(net->neurons.size());

	Stimulate(net, 0, 40);


	for (nrn_iter na = net->neurons.begin(); na != net->neurons.end(); ++na)
	{
		LANeuron* nrn = *na;
		net->spikesPerPattern[0][nrn->nid] = nrn->total_spikes;
	}

	Stimulate(net, 1, 40);


	for (nrn_iter na = net->neurons.begin(); na != net->neurons.end(); ++na)
	{
		LANeuron* nrn = *na;
		net->spikesPerPattern[1][nrn->nid] = nrn->total_spikes;
	}
	*/


	{
		int overlap=0;
		int sizeA=0;
		int sizeB=0;
		int thresh = 40;
		int activated = 0;

		for (nrn_iter ni = net->pyr_list.begin(); ni != net->pyr_list.end(); ni++)
		{
			LANeuron* n = *ni;

			if (  net->spikesPerPattern[0][n->nid] > thresh  && net->spikesPerPattern[1][n->nid] > thresh  )
				overlap++;

			if (  net->spikesPerPattern[0][n->nid] > thresh  || net->spikesPerPattern[1][n->nid] > thresh  )
				activated++;

			if ( net->spikesPerPattern[0][n->nid] > thresh) sizeA++;
			if ( net->spikesPerPattern[1][n->nid] > thresh) sizeB++;

		}

		float npyr = int( net->pyr_list.size());

		printf("Overlapping[ %f  %f ] =  %d/%d (%f %%)\n\n", 100.0*float(sizeA)/npyr, 
			100.0*float(sizeB)/npyr
			, overlap, activated ,100*float(overlap)/float(activated+0.0001) );

	}


}



// The main protocol for storing / homeostasis/recalling memories sequentially
void RunRamp(LANetwork* net)
{




	{

		ofstream rdata("./data/ramp_data.txt");
		ofstream trace("./data/ramp_voltage.txt");
		for (int i=0; i< 50; i++)
		{


			net->enablePlasticity = false;

			cout<< "ramp "<< i<<endl;

			net->injectedCurrent = i*.6;


			for (nrn_iter ni = net->pyr_list.begin(); ni != net->pyr_list.end(); ni++)
			{
				LANeuron* nrn = *ni;
				if (nrn->nid < 200)
					nrn->crebLevel =1.;
				else 
					nrn->crebLevel =0.;
			}


			net->ResetSpikeCounters();
			if (i == 40) net->traceFile = &trace;
			else net->traceFile = NULL;
			net->StimDynamics(1000);


			for (nrn_iter na = net->pyr_list.begin(); na != net->pyr_list.end(); ++na)
			{
				LANeuron* nrn = *na;
				rdata <<  nrn->total_spikes << " " ;
				cout <<nrn->total_spikes<<",";
			}
			rdata << endl;
			cout<<endl;
			
		}
	}

}






// Parse command line parameters and run model
int main( int argc, char* argv[])
{
	
	int c;


	// Set defaults
	int nneurons = 500;
	int nbranches = 20;
	int ninputs = 400;
	int nperinput = 1;
	int npatterns = 2;
	int nonesperpattern = 1;
	int interstim = 60;
	int rseed = 1980;
	char* suffix = NULL;
	bool storeData = false;
	bool disableCreb = false;

	LANetwork net; 
	net.enablePruning = false; // default
	net.runProtocol = 'S'; // Store

	while ((c = getopt(argc, argv, "M:N:H:B:I:i:P:p:T:S:s:d:w:O:g:l:h:b:c:o:t:xnLDRJCGUV"))!= -1)
	{
		switch (c)
		{
			case '?':
			cout << "usage: "<< argv[0] << " -N nneurons -B nbranches   -P npatterns -p ones_per_pattern -d overlappingPatternOffset -T interstim -S random_seed -w weakMemId " << endl;
			return 1;
			break;

			case 'B': nbranches = atoi(optarg); break;
			case 'N': nneurons = atoi(optarg); break;
			case 'P': npatterns = atoi(optarg); break;
			case 'p': nonesperpattern = atoi(optarg); break;
			case 'T': interstim = atoi(optarg); break;
			case 'S': rseed = ( atoi(optarg)); break;
			case 's': suffix = strdup(optarg); break;

			case 'x': storeData = true; break;
			case 'n': disableCreb = true; break;
			case 'w': net.isWeakMem.push_back(atoi(optarg)-1); break;

			case 'L': net.localProteins = true; break;
			case 'G': net.globalProteins = true; break;
			case 'D': net.debugMode = true; break;
			case 'R': net.repeatedLearning = true; break;
			case 'J': net.pretraining = true; break;
			case 'C': net.altConnectivity = true; break;
			case 'O': net.branchOverlap = atof(optarg); break;
			case 'H': net.homeostasisTime = atof(optarg); break;
			case 'V': net.runProtocol = 'V'; break;


			case 'o': 
				char* o = strstr(optarg, "=");
				if (o)
				{
					*o = '\0';
					char* val = o+1;

					if (!strcmp(optarg, "connectivityParam")) net.connectivityParam = atof(val); 
					else if (!strcmp(optarg,  "BSPTimeParam")) net.BSPTimeParam = atof(val); 
					else if (!strcmp(optarg,  "homeostasisTimeParam")) net.homeostasisTimeParam = atof(val); 
					else if (!strcmp(optarg,  "CREBTimeParam")) net.CREBTimeParam = atof(val); 
					else if (!strcmp(optarg,  "inhibitionParam")) net.inhibitionParam = atof(val); 

					else if (!strcmp(optarg,  "globalPRPThresh")) net.globalPRPThresh = atof(val); 
					else if (!strcmp(optarg,  "localPRPThresh")) net.localPRPThresh = atof(val); 
					else if (!strcmp(optarg,  "dendSpikeThresh")) net.dendSpikeThresh = atof(val); 

					else if (!strcmp(optarg,  "initWeight")) net.initWeight*= atof(val); 
					else if (!strcmp(optarg,  "maxWeight")) net.maxWeight*= atof(val); 
					else if (!strcmp(optarg,  "stimDurationParam")) net.stimDurationParam = atof(val); 
					else if (!strcmp(optarg,  "nNeuronsParam")) nneurons *= atof(val); 
					else if (!strcmp(optarg,  "nBranchesParam")) nbranches *= atof(val); 
					else if (!strcmp(optarg,  "nBranchesTurnover")) net.nBranchesTurnover = atoi(val); 
					else if (!strcmp(optarg,  "INClustered")) net.INClustered = atoi(val); 
					else if (!strcmp(optarg,  "setNlTypes")) net.setNlTypes = atoi(val); 
					else if (!strcmp(optarg,  "enableTurnover")) net.enableTurnover = atoi(val); 

					else if (!strcmp(optarg,  "blockedLCpattern")) net.blockedLCpattern = atoi(val); 

					else if (!strcmp(optarg,  "inSomaTau"))  net.inSomaTau  = atof(val);
					else if (!strcmp(optarg,  "pyrSomaTau")) net.pyrSomaTau = atof(val);
					else if (!strcmp(optarg,  "inDendrites")) net.inDendrites = atoi(val);

					else if (!strcmp(optarg,  "nPyr2PV")) net.nPyr2PV = atoi(val);
					else if (!strcmp(optarg,  "nPyr2SOM")) net.nPyr2SOM = atoi(val);
					else if (!strcmp(optarg,  "nPV2Pyr")) net.nPV2Pyr = atoi(val);
					else if (!strcmp(optarg,  "nSOM2Pyr")) net.nSOM2Pyr = atoi(val);

					else if (!strcmp(optarg,  "nDA2Pyr")) net.nDA2Pyr = atoi(val);
					else if (!strcmp(optarg,  "nDA2PV")) net.nDA2PV = atoi(val);
					else if (!strcmp(optarg,  "nDA2SOM")) net.nDA2SOM = atoi(val);

					else if (!strcmp(optarg,  "Pyr2SOMplastic")) net.Pyr2SOMplastic = atoi(val);
					else if (!strcmp(optarg,  "Pyr2PVplastic")) net.Pyr2PVplastic = atoi(val);

					else { 
						printf("[Error] Invalid option '%s'\n", optarg);
						exit(1);
					}

					printf("[Option] '%s' value='%s'\n", optarg, val);
				}
			break;
		}
	}


	LANetwork::SetRandomSeed(rseed);
	net.disableCreb = disableCreb;

	ninputs = nonesperpattern * npatterns;

	printf("\nInputs=%d, neurons per input=%d, memories to encode=%d\n", ninputs, nperinput, npatterns);
	net.CreateFearNet(nneurons, nbranches, ninputs, nperinput);

	char buf[1024];
	if (suffix)
		sprintf(buf, "./data/%s", suffix );
	else
		sprintf(buf, "%s", "out");

	cout << "Output data dir: "<< buf <<  endl;

	net.SetDataDir( buf);


	if (net.runProtocol == 'V')
	{
		RunRamp(&net);
	}
	else
	{
		char buf2[2048];

		interstim = 5*60;

		RunStoreTest(&net, npatterns, nonesperpattern, interstim);
		cout << "Storing data files ..."<< endl;
		cout<<buf << endl;
		net.StoreDataFiles( storeData);

		sprintf(buf2, "cp constructs.cpp %s/", buf);
		system(buf2);
		printf("Done!\n");
	}
	return 0;
}





