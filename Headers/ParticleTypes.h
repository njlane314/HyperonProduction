#ifndef _ParticleTypes_h_
#define _ParticleTypes_h_

namespace hyperon {

    inline bool isNucleus(int pdg){ return abs(pdg) > 10000; }

    inline bool isHyperon(int pdg){ return abs(pdg) == 3122 || abs(pdg) == 3212 || abs(pdg) == 3222 || abs(pdg) == 3112; }

    inline bool isNeutralSigma(int pdg){ return abs(pdg) == 3212; }
    inline bool isChargedSigma(int pdg){ return abs(pdg) == 3222 || abs(pdg) == 3112; }
    inline bool isLambda(int pdg){ return abs(pdg) == 3122; } 

    inline bool isNeutralPion(int pdg){ return abs(pdg) == 111; }
    inline bool isChargedPion(int pdg){ return abs(pdg) == 211; }
    inline bool isPion(int pdg){ return isNeutralPion(pdg) || isChargedPion(pdg); }

    inline bool isPionPlus(int pdg){ return pdg == +211; }
    inline bool isPionMinus(int pdg){ return pdg == -211; }

    inline bool isMuon(int pdg){ return abs(pdg) == 13; }

    inline bool isNucleon(int pdg){ return abs(pdg) == 2112 || abs(pdg) == 2212; }
    inline bool isProton(int pdg){ return abs(pdg) == 2212; }
    inline bool isNeutron(int pdg){ return abs(pdg) == 2112; }

    inline bool isLepton(int pdg){ return abs(pdg) == 11 || abs(pdg) == 13 || abs(pdg) == 15; }

    inline bool isNeutrino(int pdg){ return abs(pdg) == 12 || abs(pdg) == 14 || abs(pdg) == 16; }

    inline bool isChargedKaon(int pdg){ return abs(pdg) == 321; }

    inline bool isNeutralKaon(int pdg){ return abs(pdg) == 311; }
    inline bool isKaonShort(int pdg){ return abs(pdg) == 310; }
    inline bool isKaonLong(int pdg){ return abs(pdg) == 130; }

}

#endif
