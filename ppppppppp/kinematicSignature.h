#ifndef KINEMATICSIGNATURE__
#define KINEMATICSIGNATURE__
#include<string>
#include<vector>
#include"types.h"
class kinematicSignature {
	public:
		kinematicSignature () : kinematicSignature(0) {};
		kinematicSignature (size_t identifier);
		
		bool                              operator==                   (const kinematicSignature& other) const;
		realMatrix getBoseSymmetrizedKinematics (const realVector& kin)  const;
		size_t                            nKin                         ()                                const;
		sizeVector               isobarMassIndices            ()                                const;
		void                              print                        ()                                const;		
	protected:
		size_t _identifier;

};
#endif// KINEMATICSIGNATURE__
