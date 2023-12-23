/* 
    Version: $Id: constructs.cpp 172 2014-02-12 10:06:07Z gk $
    LAModel main imlementation file

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
//#include "wxglmodel.h"
#include <iterator>
#include <assert.h>


/// PROTEIN SYNTHESIS thresholds. When calcium crosses this value, proteins will be synthesized
const float GPROD_CUTOFF = 40.0; // Global ca threshold
//const float BPROD_CUTOFF = 1.8;  // Dendrite ca threshold


int LANetwork::RSEED = 123;




#define VEC_REMOVE(vec, t) (vec).erase(std::remove((vec).begin(), (vec).end(), t), (vec).end())

inline float sigmoid(float x, float x0, float broad)
{
	return (1.0 / (1.0 + exp(-x - x0)/broad));
}




// The alpha function for protein sythesis over time (x is time)

inline double nuclearproteinalpha(float x)
{
	return (x>20.)*((x-20.*60.)/(30.*60)) * exp(1. - (x-20.*60. )/(30.*60.));

	//double d =  (x-60.*60)/(40.0*60.);
	//return (exp(-d*d));

	//double d =  ((x-3.*60)/(23.*60)) * exp(1. - (x-(3.*60.) )/(23.*60.)); // 4. or 6.
	//if (d >0.) return d;
	//else return 0.;
	
}




// The alpha function for protein sythesis in dend branches over time (x is time)
inline double branchproteinalpha(float x)
{
	return ((x)/(15.*60)) * exp(1. - (x )/(15.*60.));

	
}


inline void adjust_with_bounds(float& val, float direction, float max, float min)
{
	if (direction>=0.0) val += direction*(max-val);
	else val += direction*(val - min);
}


inline float step_with_bounds(float cur_val, float dir, float max, float min)
{
	if (dir >0) return dir*(max - cur_val);
	else  return dir*(cur_val - min);
}




// the curve for the magnitude of LTP vs Ca++  . (x is calcium)
inline float caDP(float x)
{
	//return x >0.2 ? 1.0 : 0.0;
	//float f =  (2.0/(1.+exp(-(x*10.-3.5)*10.))) - (1.0/(1.+exp(-(x*10.0-0.5)*19.)));


	float f = (1.3/(1.+exp(-(x*10.-3.5)*10.))) - (0.3/(1.+exp(-(x*10.0-2.0)*19.)));
	return f;
	//return f;//(2.0/(1.+exp(-(x*10.-0.7)*10.))) - (1.0/(1.+exp(-(x*10.0-0.2)*19.)));
	//return (2.0/(1.+exp(-(x*10.-3.1)*8.))) - (1.0/(1.+exp(-(x*10.0-2.1)*6.)));
}







// Create a list of neurons
void LANetwork::CreateNeurons(int number, int n_branches_per_neuron, char type, vector<LANeuron*>* appendTo, int inputId, int somethingDummy)
{
	for (int i =0 ;i < number; i++)
	{
		LANeuron* n ;
		if (type == 'S')
		{
			n = new LAInput();
			n->network = this;

			n->input_id = inputId;
			((LAInput*)n)->groupIdx = i;
		}
		else
		{
			n = new LANeuron();
			n->network = this;
			for (int tt =0; tt < n_branches_per_neuron; tt++)
			{
				LABranch* bb  = new LABranch;
				bb->bid = this->branches.size();
				bb->neuron = n;

				if (type == 'P')  // Pyramidals
				{
					bb->nlType = DEND_SUPRA;
				}
				else if ( type == 'M') // SOM interneurons
				{
					bb->nlType = DEND_LINEAR;

				}
				else if (type == 'V') // basket interneurons
				{
					bb->nlType = DEND_LINEAR;
				}


				this->branches.push_back(bb);
				n->branches.push_back(bb);
			}
		}
		n->type = type;

		n->V = param_V_rest +1.0;
		n->nid = this->neurons.size();

		n->glx = 1.9*(n->nid%50)/50. - 0.9;
		n->gly = 1.9*(n->nid/50)/50. - 0.9;

		this->neurons.push_back(n);
		if (appendTo)
			appendTo->push_back(n);
	}
}




void LANetwork::CalculateDistances()
{
	for (nrn_iter na = neurons.begin(); na != neurons.end(); ++na)
		for (nrn_iter nb = neurons.begin(); nb != neurons.end(); ++nb)
			if (nb != na)
			{
				LANeuron* a = *na, *b=*nb;

				
				float dx = b->pos_x - a->pos_x;
				float dy = b->pos_y - a->pos_y;
				float dz = b->pos_z - a->pos_z;
				distances[pair<int,int>(a->nid,b->nid)] = (double)sqrt(dx*dx + dy*dy + dz*dz);
			}
}


/*
int LANetwork::ConnectToRandomBranches(vector<LANeuron*> fromList, vector<LANeuron*> toList,   int nSynapses, float weight)
{
		LANeuron* from = fromList.at(int(rgen()*(float)toList.size()));
		LANeuron* to = toList.at(int(rgen()*(float)toList.size()));

		LASynapse* syn = new LASynapse();
		syn->sid  = this->synapsesCounter++;
		syn->target_branch = b;

		syn->source_nrn = from;
		syn->target_nrn = to;
		syn->isPlastic = isPlastic;
		syn->weight  = weight;

		syn->target_branch = syn->target_nrn->branches[(int)rval];

		syn->source_nrn->outgoing.push_back(syn);
		syn->target_nrn->incoming.push_back(syn);
		syn->target_branch->synapses.push_back(syn);
		syn->pos = rgen();
		this->synapses.push_back(syn);

}

*/


void LANetwork::AddSynapse(LANeuron* a, LABranch* br, float weight, bool isPlastic)
{
		LASynapse* syn = new LASynapse();
		syn->sid  = this->synapsesCounter++;
		syn->source_nrn = a;
		syn->target_nrn = br->neuron;
		syn->isPlastic = isPlastic;

		syn->weight  = weight;
		syn->target_branch = br; 

		syn->source_nrn->outgoing.push_back(syn);
		syn->target_nrn->incoming.push_back(syn);
		syn->target_branch->synapses.push_back(syn);
		syn->pos = rgen();
		this->synapses.push_back(syn);
}



int LANetwork::ConnectNeurons(vector<LANeuron*> fromList, vector<LANeuron*> toList, bool isClustered, float toDistance, int nNeuronPairs, int nSynapsesPerNeuron, float weight, bool isPlastic , bool randomizeweight , float overlap)
{
	int tpairs =0;
	if (nNeuronPairs <=0) return 0;
	while(true)
	{
		LANeuron* a = fromList.at(int(rgen()*(float)fromList.size()));
		LANeuron* b = toList.at(int(rgen()*(float)toList.size()));

		for (int i =0; i < nSynapsesPerNeuron; i++)
		{

			float rval;
			if (isClustered)
				rval = rgen()*float(b->branches.size()/3); /// XXX Clustering
			else
				rval = rgen()*float(b->branches.size());

			LABranch* br =  b->branches[(int)rval];

			this->AddSynapse(a, br, weight, isPlastic);	
		}

		if (++tpairs >= nNeuronPairs) break;
	}
	return tpairs;
}




inline int randomindex(int max)
{
	return (int)(rgen()*float(max));
}

int LANetwork::PurgeInputSynapses(int totalToRemove,float weightLimit)
{

	int totalFound =0;
	int totalTries = totalToRemove*3;
	while (totalFound < totalToRemove && totalTries-->0)
	{
		nrn_list lst = this->inputs_cs.at(randomindex(this->inputs_cs.size()));
		LANeuron* n = lst.at(randomindex(lst.size()));
		if (n->outgoing.size())
		{
			LASynapse* s = n->outgoing.at(randomindex(n->outgoing.size()));
			if (s->target_nrn->type == 'P'  && s->weight <= weightLimit)
			{
				//candidate for deletion

				//pthread_mutex_lock(&this->synapses_mutex);

				VEC_REMOVE(s->source_nrn->outgoing, s);
				VEC_REMOVE(s->target_nrn->incoming, s);
				VEC_REMOVE(s->target_branch->synapses, s);
				VEC_REMOVE(this->synapses, s);
				//cout << " Removing " << s->sid << endl;
				delete s; 

				//pthread_mutex_unlock(&this->synapses_mutex);

				totalFound++;
			}
		}
	}

	if (totalFound < totalToRemove)
	{
		cout << " Warning: not enougn synapses to remove: " << totalToRemove << endl;  
	}

	return totalFound;
}




void LANetwork::CreateFearNet(int nneurons, int nbranches, int ninputs, int nneuronsperinput)
{
	this->n_neurons = nneurons;
	this->n_inputs = ninputs;
	this->n_neurons_per_input = nneuronsperinput;
	this->n_branches_per_neuron = nbranches;

	// Create Pyrs
	this->CreateNeurons(this->n_neurons*0.8, this->n_branches_per_neuron, 'P', &this->pyr_list, -1, this->nBranchesTurnover);

	this->CreateNeurons(this->n_neurons*0.1, this->inDendrites , 'V', &this->in_pv); // PV

	this->CreateNeurons(this->n_neurons*0.1, this->inDendrites , 'M', &this->in_som); // SOM

	
	for (int i=0;i < n_inputs;i++)
	{
		vector<LANeuron*> nrnsCS;
		CreateNeurons(this->n_neurons_per_input, 0, 'S', &nrnsCS, i, true);
		this->inputs_cs.push_back(nrnsCS);
	}


	CreateNeurons(10, 0, 'S', &this->da_list, -1, true);

	this->CalculateDistances();

	float synBase = 1.0;
	//float synW = 1.;


	//xcon connectivity
	if (runProtocol != 'V')
	{
	ConnectNeurons(this->pyr_list, this->in_pv,    0, 10.,  synBase*  1000, 1, .6, false); 
	ConnectNeurons(this->in_pv,    this->pyr_list, 0, 10.,  synBase* 10000, 1, 0.5, false); 

	ConnectNeurons(this->pyr_list, this->in_som, 0, 10., synBase* 1000, 1, .3, false); 
	ConnectNeurons(this->in_som, this->pyr_list, 0, 10., synBase*2000, 1, .3, false); 
	}


	for (int i =0; i < this->n_inputs ; i++)
	{
		this->ConnectNeurons(this->inputs_cs[i], this->pyr_list, 0, 10.0, 5600, 1, .4 /* initWeight */, true, false, this->branchOverlap);
		//this->ConnectNeurons(this->inputs_cs[i], this->in_pv, 0, 10.0,1000, 1, .5 /* initWeight */, false, false, this->branchOverlap);
		//this->ConnectNeurons(this->inputs_cs[i], this->in_som, 0, 10.0,1000, 1, .5 /* initWeight */, false, false, this->branchOverlap);
	}


	// Dopamine stimulation
	if (this->nDA2Pyr>100)
	{
		ConnectNeurons(this->da_list, this->pyr_list,  0.0, 10., 400, 1, 1., false);
		//ConnectNeurons(this->da_list, this->in_pv,     0.0, 10., 8000, 1, 1., false);
		//ConnectNeurons(this->da_list, this->in_som,    0.0, 10., 8000, 1, 1., false);
	}

	this->RecordInitialWeightSums();
}











template <typename T> static void PrintVector( vector<T>&  ar, ostream& outfile) 
{
	for (typename vector<T>::iterator it = ar.begin(); it != ar.end(); it++)
	{
		outfile << *it << ' ';
	}
	outfile << std::endl;
}





void LANetwork::StoreDataFiles( bool extras = false )
{
	vector<int> pattern;

	string dirname = this->datadir;
	ofstream paramsdat((dirname + "/parameters.txt").c_str());
	paramsdat <<"total_neurons="<< this->neurons.size() << endl;
	paramsdat <<"total_pyramidals="<< this->pyr_list.size() << endl;
	paramsdat <<"branches_per_neuron="<< this->n_branches_per_neuron << endl;
	paramsdat <<"number_inputs="<< this->inputs_cs.size() << endl;
	paramsdat <<"neurons_per_input="<< this->n_neurons_per_input << endl;
	paramsdat <<"rseed="<< RSEED << endl;


	ofstream patternsdat((dirname + "/patterns.txt").c_str());
	for (vector<vector<int> >::iterator it = this->patterns.begin(); it != this->patterns.end(); it++)
	{
		pattern = *it;
		copy(pattern.begin(), pattern.end(), ostream_iterator<int>(patternsdat, " "));
		patternsdat << endl;
	}

	ofstream synstatedat((dirname + "/synstate.dat").c_str());

	ofstream spikesdat((dirname + "/spikes.dat").c_str());
	ofstream crebdat((dirname + "/creb.dat").c_str());
	ofstream voltagedat((dirname + "/voltages.dat").c_str());
	ofstream branchspikesdat((dirname +"/branchspikes.dat").c_str());
	ofstream branchcalcium((dirname + "/branchcalcium.dat").c_str());
	ofstream weightsdat((dirname + "/weights.dat").c_str());
	ofstream branchproteins((dirname + "/branchproteins.dat").c_str());
	ofstream branchstrengths((dirname + "/branchstrengths.dat").c_str());
	ofstream tagsdat((dirname + "/tags.dat").c_str());
	ofstream nrnproteindat((dirname + "/nrnprotein.dat").c_str());
	ofstream weighthistorydat((dirname + "/weighthistory.dat").c_str());
	ofstream dbgneuron((dirname + "/dbgneuron.dat").c_str());

	ofstream syn_per_branch((dirname + "/syn_per_branch.dat").c_str());

	PrintVector<float>( dbgNeuron, dbgneuron);

	for (nrn_iter na = neurons.begin(); na != neurons.end(); ++na)
	{
		LANeuron* nrn = *na;

		PrintVector<int>( nrn->spikeTimings, spikesdat);
		PrintVector<float>( nrn->proteinHistory, nrnproteindat);

		for (branch_iter bi = nrn->branches.begin(); bi != nrn->branches.end(); ++bi)
		{
			LABranch* b = *bi;

			PrintVector<float>(b->branchSpikesHistory, branchspikesdat);

			for (syn_iter si = b->synapses.begin(); si != b->synapses.end(); ++si)
			{
				LASynapse* s = *si;
				//PrintVector<float>(s->weightHistory, weightsdat);
				synstatedat << s->sid<<" " 
					<< b->bid<<" "
					<< nrn->nid  << " "
					<< s->source_nrn->nid <<" " 
					<< s->source_nrn->input_id<< " "
					<< b->strength  << " "
					<< s->weight << " " 
					<<endl;

				if (s->tagHistory.size())
				{
					tagsdat << s->source_nrn->input_id<< " ";
					//PrintVector<float>(s->tagHistory, tagsdat);
				}

				if (s->weightHistory.size())
				{
					weighthistorydat << s->source_nrn->input_id<< " ";
					//PrintVector<float>(s->weightHistory, weighthistorydat);
				}

				if (s->isPlastic && s->source_nrn->input_id >=0 && s->weight > .7)
				{
				//	totPot[ s->source_nrn->input_id >=0 ] += 1;
				}
			}

		}
	}

	ofstream sppdat((dirname + "/spikesperpattern.dat").c_str());
	for (uint i =0; i < this->spikesPerPattern.size(); i++)
	{
		for (uint j =0; j < this->spikesPerPattern[i].size(); j++) sppdat << this->spikesPerPattern[i][j] << " ";
		sppdat << endl;
	}



}




// The main loop that simulates the stimulation dynamics, voltages,  Ca influx etc. 
// The time step is 1msec 
//
void LANetwork::StimDynamics(int duration)  // duration in msec
{
	int t = 0;
	//bool spikeState[this->neurons.size()+1];
	//int lastSpikeT[this->neurons.size()+1];
	//int totSpikes=0; 
	


	// Zero out spike states
	//fill_n(spikeState, this->neurons.size()+1, 0);
	//fill_n(lastSpikeT, this->neurons.size()+1, 0);

	// zero out calcium
	for (syn_iter si = this->synapses.begin(); si != this->synapses.end(); ++si)
	{
		(*si)->calcium = 0.0;
	}

	// zero out neuron voltages
	for (nrn_iter ni = this->neurons.begin(); ni != this->neurons.end(); ++ni)
	{
		LANeuron* n = *ni;
		n->wadapt = 0.0;
		n->vreset =0.;
		n->bAP =0;
		n->lastSpikeT = -100;
		//n->iCaap =0.;
		//n->Vfilt =0.;
		n->V =0.;
		n->totcalc =0.;
	}


	// zero out branch calcium counters + voltages
	for (branch_iter bi=this->branches.begin(); bi != this->branches.end(); ++bi)
	{
		(*bi)->totcalc = 0.0;
		(*bi)->depol = 0.0;
		(*bi)->dspike = 0.0;
		(*bi)->dspikeT = -1;
	}

	const float g_ampa = 0.13;
	const float g_gaba = 0.1;

	float soma_exc =0;  // soma  current
	float soma_gaba =0;  // soma inhibition


	// Run simulation steps for each msec


	for (t=0; t < duration; t++)
	{
		for (nrn_iter ni = this->neurons.begin(); ni != this->neurons.end(); ++ni)
		{
			LANeuron* n = *ni;

			if (n->type == 'S') // Source/stimulation neuron
			{
				LAInput* in = (LAInput*)n;
				if (in->spikeTimes && t >= in->nextSpikeT && in->nextSpikeT >0)
				{
					// Time to emit a spike
					if (in->curSpike < in->totalSpikes)
						in->nextSpikeT = in->spikeTimes[in->curSpike++];
					else 
						in->nextSpikeT = -1;

					//spikeState[in->nid] = 1;
					//n->isSpiking = true;
					//lastSpikeT[in->nid] = t;
					n->lastSpikeT = t;
					in->total_spikes++;
					in->V = 20.0; 
				}
				else
				{
					//spikeState[in->nid] = 0;
					//n->isSpiking = false;
					//in->V = 0.0;
				}
			}
			else
			{

				soma_exc =0;  // soma  current
				soma_gaba =0;  // soma inhibition

				// dive into all branches and calculate local depolarization / calcium etc
				for (branch_iter bi=n->branches.begin(); bi != n->branches.end(); ++bi)
				{
					LABranch* b = *bi;
					float dend_exc =0.;
					float dend_inh =0.;

					if (n->lastSpikeT == t-1) b->depol =0;

					// gather input spikes
					for (syn_iter si=b->synapses.begin(); si != b->synapses.end(); ++si)
					{
						LASynapse* s = *si;
						if (s->source_nrn->lastSpikeT == t-1
						 || s->source_nrn->lastSpikeT == t-2)
						{
							// PV synapses are connected to branches, but in reality they are to be
							// redirected to the soma
							if (s->source_nrn->type == 'V' ) soma_gaba += ( s->weight);
							else if (s->source_nrn->type == 'M' ) dend_inh += ( s->weight);
							else dend_exc += (s->weight);
						}
					}


					if (b->nlType == DEND_LINEAR)
					{
						// linear  interneuron dendrites 
						b->depol = b->depol + dend_exc*g_ampa*(70-b->depol)  - dend_inh*g_gaba*(10+b->depol) -  b->depol/20.0;
					}
					else  
					{
						// excitatory neuron
						b->depol = b->depol + dend_exc*g_ampa*(70-b->depol) 
							- dend_inh*g_gaba*(10+b->depol) 
							- b->depol/25.0;
						
					}


					if (b->depol > 70) b->depol = 70;
					else if(b->depol < -10) b->depol = -10;

					//   calcium influx from synapses
					if (this->enablePlasticity && t>0)
					{
						for (syn_iter si=b->synapses.begin(); si != b->synapses.end(); ++si)
						{
							LASynapse* s = *si;
							if (s->isPlastic && s->source_nrn->lastSpikeT == t-1)
							{
								if (n->bAP>1) //dep >30)
								{
									float ff = .50; 
									s->calcium += ff;
								}
							}
						}
					}

					float df = b->depol - n->V ;
					if (df>0)
						soma_exc += ( 0.15 )*df;

				}

				if ( t < n->lastSpikeT + 3 ) // Refractory period
				{
					if ( t == n-> lastSpikeT + 2 ) // End of refractoriness
					{
						n->V =-3; // Reset
					}
					else if (t == n->lastSpikeT +1) // Just emitted spike
					{
						n->V = -3;
						if (n->crebLevel>0) n->wadapt += 1.0+ n->wadapt*1.0; //*1.4 ; //+ 4.0;
						else 
							n->wadapt += 2.8 +  n->wadapt*.45; //*1.4; //+ 9.0*2;
					}
				}
				else
				{
					if (n->type == 'V' || n->type == 'M') //  interneuron
					{
						n->V +=  soma_exc -  (n->V)/20.; 
					}
					else
					{
						// Pyr neuron
						float noisecur =  rgen()*(0.5  +this->injectedCurrent) ;
						n->somaInh += soma_gaba - n->somaInh/30;

						n->V +=  noisecur + soma_exc - n->somaInh*3.0 - (n->V)/20. - n->wadapt;  

						// CREB activation changes the adaptation dynamics
						float tc = 210;
						if (n->crebLevel>0)  tc  = 100.;

						n->wadapt +=  ( 0.02*n->V -  n->wadapt)/( tc); 

					}


					if (n->V > 70) n->V = 70;
					else if (n->V < -10) n->V = -10;

					float thr = 18.0;

					if (  n->V > thr)
					{
						// Generate a somatic spike
						n->lastSpikeT = t;
						n->total_spikes++;
						n->V = 70.;
						n->bAP = 5.;
						n->totcalc += 1;

					}

				}

				if (n->bAP>0) n->bAP-=n->bAP / 10; 
				// adaptation cap
				if (n->wadapt <0.) n->wadapt =0.;

			}


			if (n->lastSpikeT == t)
			{
				// Record this spike
				n->spikeTimings.push_back(t+ this->Tstimulation);
			}

			if (this->traceFile != NULL)
			{
				if (n->nid == 2) 
					(*this->traceFile) << n->V << " ";
				else if (n->nid == 202)
					(*this->traceFile) << n->V << endl;
			}


		}

	}

	this->Tstimulation += duration;
}




void LANetwork::SaveSnapshot(char* filename)
{
	ofstream of(filename);
	for (nrn_iter ni = this->pyr_list.begin(); ni != this->pyr_list.end(); ni++)
	{
		LANeuron* nrn = *ni;
		if (nrn->type == 'P')
		{
			float totTag =0;
			float totProtein =0;
			for (branch_iter bi = nrn->branches.begin(); bi != nrn->branches.end(); ++bi)
			{
				LABranch* b = *bi;
				for (syn_iter si = b->synapses.begin(); si != b->synapses.end(); ++si)
				{
					LASynapse* s = *si;
					if (s->isPlastic)
					{
						totTag += s->stag;
					}
				}
				totProtein += b->protein;
			}
			of << nrn->nid << ' ' << totProtein << ' '<<  totTag << endl;
		}
	}

}



void LANetwork::Interstim(int durationSecs)
{
	int tstop = T + durationSecs;
	this->isInterstim = true;

	printf("Interstim %d seconds (T=%d) plast=%d G=%d, L=%d ... \n", durationSecs, T, this->enablePlasticity, this->globalProteins, this->localProteins);

	float tstep = 60.0;


	int trec =0, tactrec=0;
	int totTagged=0, totTaggedMinus=0, totBProd =0;
	float totLTP=0., totLTD=0.;
	int totact =0;
	int totbspikes =0, totSpikes=0;
	float maxSpk =0;
	

	for (nrn_iter ni = this->pyr_list.begin(); ni != this->pyr_list.end(); ni++)
	{
		LANeuron*  nrn = *ni;
		float nrnCalc =0.0;
		for (branch_iter bi = nrn->branches.begin(); bi != nrn->branches.end(); ++bi)
		{
			LABranch* b = *bi;
			totbspikes += b->branch_spikes;

			if (!this->enablePlasticity)
				continue;

			for (syn_iter si = b->synapses.begin(); si != b->synapses.end(); ++si)
			{
				LASynapse* s =*si;

				if (s->calcium > 1.0)
				{
					if (nrn->total_spikes > GPROD_CUTOFF)
					{
						s->stag = 1.0;
						totTagged++;
					}
					else 
					{
						s->stag = -0.3;
						totTaggedMinus++;
					}
				}


				b->totcalc += s->calcium;
				s->calcium = 0.;
			}


			
			if ( b->totcalc  > this->localPRPThresh) // This branch should produce PRPs now BPROD
			{
				b->prpTransients.push_back( pair<float,float>(T, b->totcalc));
				totBProd++;
			}

			nrnCalc +=  b->totcalc;

		}

		totSpikes += nrn->total_spikes;

		if (nrn->total_spikes > CUTOFF*4.0)
		{
			totact++;
		}

		if (maxSpk < nrn->total_spikes) maxSpk = nrn->total_spikes;


		if (this->enablePlasticity)
		{
			if (nrn->totcalc > GPROD_CUTOFF) // Protein synthesis trigger
			{
				nrn->prpTransients.push_back( pair<float,float>(T, nrnCalc) );

				if (nrn->total_spikes > CUTOFF*4)
					tactrec ++;
				trec ++;
				if (!this->disableCreb ) nrn->crebLevel=1.0;
			}
			else if (!this->disableCreb ) nrn->crebLevel=-1.0;
		}

	}

	
	
	printf("Tagged: [%d/%d +%.1f/%.1f] Gprod:%d (%.1f%%) Bprod:%d, Act:%d (%.1f%%), Act&Gprod:%d Freq:%.1f (max %.1f) Dspikes:%d\n\n", totTagged, totTaggedMinus, totLTP, totLTD, trec, 100.*(float)trec/(float(this->pyr_list.size())), totBProd, totact, 100.*(float)totact/(float(this->pyr_list.size())), tactrec, (float)totSpikes/((float)this->pyr_list.size()*4.0), float(maxSpk)/4.0, totbspikes );

	for (; T < tstop; T+= tstep)
	{
		for (nrn_iter ni = this->pyr_list.begin(); ni != this->pyr_list.end(); ni++)
		{
			LANeuron*  nrn = *ni;

			float totalSynapseWeight =0.0;
			float totalBranchStrength =0.0;
			int totalSynapses =0;
			nrn->protein =0.0;

			// "Whole-neuron" distribution of proteins
			nrn->proteinRate =0;
			for (pair_iter ii = nrn->prpTransients.begin(); ii != nrn->prpTransients.end(); ++ii)
			{
				pair<float, float> p = *ii;
				int td= (T - p.first);
				float al = (nuclearproteinalpha(td));
				if (nrn->proteinRate < al)
					nrn->proteinRate =  al;
			}

			
			
			for (branch_iter bi = nrn->branches.begin(); bi != nrn->branches.end(); ++bi)
			{
				LABranch* b = *bi;
				b->proteinRate =0.;

				for (pair_iter ii = b->prpTransients.begin();ii != b->prpTransients.end(); ++ii)
				{
					pair<float, float> p = *ii;
					float td = float(T - p.first);
					float al = (branchproteinalpha(td));
					if (b->proteinRate < al)
						b->proteinRate = al;
				}


				float  f =0.;

				if (this->localProteins)
					f = b->proteinRate; 
				else if (this->globalProteins)
					f =  nrn->proteinRate;
				else
				{
					f = 1.0*b->proteinRate + 1.0* nrn->proteinRate;
					if (f>1.0) f = 1.0;
				}



				b->protein = this->proteinsParam* f;

				vector<LASynapse*> candidates;


				for (syn_iter si = b->synapses.begin(); si != b->synapses.end(); ++si)
				{
					LASynapse* s =*si;

					if (s->stag != 0.0) 
					{
						s->stag -= (tstep/3600.)* s->stag;

						if (b->protein > 0.1 && (s->stag >0.1 || s->stag < 0.1))
						{
							//candidates.push_back(s);
							float fw = s->stag* b->protein;
							//if (s->stag >0.)
							s->weight += tstep * fw/400.;
						}
					}

					if(1)
					{
						if (s->weight  > maxWeight)
						{
							s->weight = maxWeight;
							//s->stag =0.;
						}
						else if (s->weight  < 0.)
						{
							s->weight = 0.;
							//s->stag =0.;
						}
					}


					totalSynapseWeight += s->weight;
					totalSynapses++;

					//Homeostasis
			
					s->weight += s->weight * (1.0 - nrn->synScaling )*tstep/(30.*24.0* 3600.*homeostasisTimeParam); // Synaptic scaling 


				}

				



				if (T%800 ==0)
				{
					b->branchProteinHistory.push_back(b->protein);
				}
			}


			// Synaptic  scaling
			if (totalSynapses>0)
				nrn->synScaling = totalSynapseWeight / (initWeight*float(totalSynapses));
		
			else
				nrn->synScaling = 1.0;

			// Branch plasticity homeostasis
			nrn->branch_scaling = totalBranchStrength/((float(1.0) * float(nrn->branches.size())));


		
		}


		#ifdef WXGLMODEL_H
		if (this->wx_window && T%20 ==0)
		{
			this->wx_window->UpdatePanel();
		}
		#endif
	}
	



	if (this->enablePlasticity && this->enableTurnover) this->DoTurnover(durationSecs);

	this->isInterstim = false;
}



void LANetwork::DoTurnover(float durationSecs )
{

	float pTurnOver = 4.*durationSecs/86400.; // per day

	int synsAdded=0;
	for (branch_iter bi = this->branches.begin(); bi != this->branches.end(); bi++)
	{
		LABranch* b = *bi;
		//if (b->turnoverRate>0.0)
		{
			vector<LASynapse*> v;
			for (syn_iter si = b->synapses.begin(); si != b->synapses.end(); si++)
			{
				LASynapse* s = *si;
				if (s->isPlastic && s->weight <= this->initWeight)
				{
					float p;
					{
						float dist = s->pos - b->turnoverPoint;
						p = exp(-(dist*dist)/.1)*b->turnoverRate;
					}

					if (rgen() < p*pTurnOver)
					{
						VEC_REMOVE(s->source_nrn->outgoing, s);
						VEC_REMOVE(s->target_nrn->incoming, s);
						VEC_REMOVE(this->synapses, s);
						//cout << " Removing syn " << s->sid << endl;

						si = b->synapses.erase(si); //VEC_REMOVE(s->target_branch->synapses, s);

						delete s; 

						/* Connect  a random input back here */

						nrn_list lst = this->inputs_cs.at(randomindex(this->inputs_cs.size()));
						LANeuron* n = lst.at(randomindex(lst.size()));
						this->AddSynapse(n, b, initWeight, true);
						synsAdded++;

						continue;
					}
					
				}
			}

		}
	}

	printf("Added/removed %d synapses\n", synsAdded);
		

}



void LANetwork::Begin()
{
	
}



void LANetwork::RunTests()
{

	LANetwork net; 

	vector<LANeuron*> pyr;
	vector<LANeuron*> inh;
	vector<LANeuron*> inp;

	net.CreateNeurons(10, 10, 'P', &pyr, -1, 0);
	cout << (pyr.size() == 10) << endl;

	net.CreateNeurons(10, 10, 'V', &inh, -1, 0);
	cout << (inh.size() == 10) << endl;

	net.CreateNeurons(10, 0, 'S', &inp, -1, 0);
	cout << (inp.size() == 10) << endl;

	for (nrn_iter na = pyr.begin(); na != pyr.end(); ++na)
	{
		LANeuron* nrn  = *na;;
		cout << (nrn->type == 'P') << endl;
	}

	cout << (net.synapses.size() ==0 ) << endl;

	net.ConnectNeurons(pyr, inh, false, (float)10.,  1000, 1, 1.0, false);

	cout << (net.synapses.size()  == 1000) << endl;
}



