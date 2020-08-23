

// Forward declarations
class Model;
class HJB;

class StationaryDist {
	public:
		StationaryDist() {}

		void compute(const Model& model, const HJB& hjb);

};