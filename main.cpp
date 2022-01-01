#include <mod/Chem.hpp>
#include <mod/graph/Graph.hpp>
#include <mod/graph/internal/Graph.hpp>
#include <mod/rule/Rule.hpp>
#include <mod/rule/internal/Rule.hpp>
#include <mod/lib/Chem/MoleculeUtil.hpp>
#include <mod/lib/Graph/Single.hpp>
#include <mod/lib/Graph/Properties/Molecule.hpp>
#include <mod/lib/Graph/LabelledGraph.hpp>
#include <mod/lib/LabelledUnionGraph.hpp>
#include <mod/lib/Rules/Real.hpp>

#include <boost/bimap.hpp>

#include <iostream>

// import asRange which wraps a std::pair of iterators with a proxy class
// to enable ADL of begin() and end(), so it can be used with ranged-based for loops
using mod::asRange;

using GraphType = mod::lib::Graph::GraphType;
using Vertex = mod::lib::Graph::Vertex;
using Edge = mod::lib::Graph::Edge;
using LUG = mod::lib::LabelledUnionGraph<mod::lib::Graph::LabelledGraph>;
using UVertex = boost::graph_traits<LUG::GraphType>::vertex_descriptor;
using VertexMap = boost::bimap<UVertex, UVertex>;
using AtomData = mod::AtomData;
using BondType = mod::BondType;
using PropString = mod::lib::Graph::LabelledGraph::PropStringType;
using MolView = mod::lib::Graph::PropMolecule;

using mod::lib::Chem::bondToChar;

std::shared_ptr<mod::rule::Rule> createRule(const VertexMap &vertexMap,
                                            const LUG &lgEduct, const LUG &lgProduct, bool doChemistryCheck);

template<typename Graph>
std::size_t getVertexId(const typename boost::graph_traits<Graph>::vertex_descriptor &v, const Graph &g) {
	return get(boost::vertex_index_t(), g, v);
}

template<typename Graph>
auto getVertexFromId(std::size_t id, const Graph &g) {
	assert(id < num_vertices(g));
	return vertex(id, g);
}

int bondValue(BondType bt) {
	switch(bt) {
	case BondType::Single:
		return 1;
	case BondType::Double:
		return 2;
	case BondType::Triple:
		return 3;
	case BondType::Aromatic:
		throw mod::InputError("Aromatic bonds are not supported.");
	}
	__builtin_unreachable();
}

std::vector<std::shared_ptr<mod::rule::Rule> > doStuff(
		const std::vector<std::shared_ptr<mod::graph::Graph>> &educts,
		const std::vector<std::shared_ptr<mod::graph::Graph>> &products,
		bool doChemistryCheck) {
	// first make objects representing the disjoint union of the educts and products
	const auto makeUnion = [](const auto &gs) {
		LUG lug;
		for(const auto &g : gs)
			mod::graph::internal::push_back(lug, &mod::graph::internal::getLabelledGraph(g->getGraph()));
		return lug;
	};
	const auto lgEduct = makeUnion(educts);
	const auto lgProduct = makeUnion(products);
	// the underlying graphs of the unions
	const auto &gEduct = get_graph(lgEduct);
	const auto &gProduct = get_graph(lgProduct);

	if(num_vertices(gEduct) != num_vertices(gProduct)) {
		throw mod::InputError("Error: the graphs do not have the same number of vertices:\n\teduct: "
		                      + std::to_string(num_vertices(gEduct)) + "\tproduct: "
		                      + std::to_string(num_vertices(gProduct)) + "\n");
	}

	{ // examples of how to access stuff	
		std::cout << "Educt graph:\n";
		for(const auto v : asRange(vertices(gEduct))) { // iterating through vertices of a graph
			// Getting the ID of a vertex (in the range [0; num_vertices[).
			// getVertexId is defined as a short-hand in this file.
			std::cout << getVertexId(v, gEduct)
			          // Getting the associated data on a vertex, when it is interpreted as an atom.
			          // mod::graph::internal::getMolecule() on a vertex returns a mod::AtomData
			          // See: http://jakobandersen.github.io/mod/libmod/Chem.html#_CPPv4N3mod8AtomDataE
			          << "\t\"" << mod::graph::internal::getMolecule(v, lgEduct) << "\"\n";
			for(const auto e : asRange(out_edges(v, gEduct))) { // iterating through the incident edges of a vertex
				// Representing an undirected edge as a set of two vertices would provide a bad interface.
				// So the endpoints are distinguished as the 'source' and the 'target' of the edge,
				// even though the edge is undirected. There are thus two edge descriptors that represent the same
				// undirected edge, with the source and target swapped. Such two descriptors still compare equal with operator==.
				// Therefore it makes sense to not have a single 'incident_edges()' function, but instead two different
				// functions out_edges() and in_edges() which both iterate through the incident edges, but with source always
				// being the query vertex for out_edges() and target for in_edges().
				const auto t = target(e, gEduct); // so, here we get the other endpoint than v
				// Access the associated data on an edge, interpreted as if it models a bond in a molecule.
				// mod::graph::internal::getMolecule() on an edge returns a mod::BondType
				// See: http://jakobandersen.github.io/mod/libmod/Chem.html#_CPPv4N3mod8BondTypeE
				std::cout << "\t\"" << mod::graph::internal::getMolecule(e, lgEduct)
				          << "\"\t" << getVertexId(t, gEduct) << '\n';
			}
		}

		std::cout << "Product graph:\n";
		// ... and all this is the exact same as for the educt graph, just with the product graph
		for(const auto v : asRange(vertices(gProduct))) {
			std::cout << getVertexId(v, gProduct) << "\t\"" << mod::graph::internal::getMolecule(v, lgProduct) << "\"\n";
			for(const auto e : asRange(out_edges(v, gProduct))) {
				const auto t = target(e, gProduct);
				std::cout << "\t\"" << mod::graph::internal::getMolecule(e, lgProduct)
				          << "\"\t" << getVertexId(t, gProduct) << "\n";
			}
		}
	}

	std::vector<VertexMap> vertexMaps;
	{ // START WORKING HERE
		{ // Example ONE --- hard-coding maps
			// This probably only works when the example graphs are loaded, exactly as they are written.

			// Make a mapping based on accessing the vertices via their IDs - not a good idea!
			// Let's nevertheless do it to have a starting point to improve:

			// First a little helper-function that converts vertex IDs into the corresponding descriptors,
			// and then adds to the vertex map.
			const auto setById = [&gEduct, &gProduct](VertexMap &vm, std::size_t idEduct, std::size_t idProduct) {
				vm.insert(VertexMap::value_type(getVertexFromId(idEduct, gEduct), getVertexFromId(idProduct, gProduct)));
			};

			{ // first map
				VertexMap vm;
				for(std::size_t i = 0; i < 6; i++)
					setById(vm, i, i);
				setById(vm, 6, 7);
				setById(vm, 7, 6);
				vertexMaps.push_back(std::move(vm));
			}


			{ // second map
				VertexMap vm;
				setById(vm, 1, 1);
				setById(vm, 2, 2);
				setById(vm, 3, 3);
				setById(vm, 6, 7);
				vertexMaps.push_back(std::move(vm));
			}

			{ // third map, intended to be a duplicate of the second map
				VertexMap vm;
				setById(vm, 1, 1);
				setById(vm, 2, 2);
				setById(vm, 3, 3);
				setById(vm, 6, 7);
				vertexMaps.push_back(std::move(vm));
			}
		} // end of Example ONE

		{ // Example TWO --- an elaboration of Example ONE
			VertexMap vm;
			// This is just fo explanation of what setById is really doing,
			// using more elaborate code.
			//
			// The following code illustrates how, as an example,
			//
			// setById(vm, 7, 6);
			//
			// could be implemented using our own implementation of getVertexFromId() as a loop,
			// to find the proper vertex descriptors for the vertices with IDs 6 and 6 in their respective graphs.
			// The variables "seven" and "six" are those descriptors.
			//
			// As vertex IDs are assigned by how a graph was originally constructed, it is not a good idea to use them
			// for anything else than indices into arrays (or debugging).
			// As the code should not be executed here it is wrapped in "if(false)", but you
			// can literally replace "setById(vm, 7, 6)" from the example above by the code below.

			if(false) {
				// --- REPLACE setById(vm, 7, 6); START ---
				UVertex seven, six;
				for(const auto v : asRange(vertices(gEduct))) // "for each vertex v in gEduct:
					if(getVertexId(v, gEduct) == 7)
						seven = v;
				for(const auto v : asRange(vertices(gProduct))) // "for each vertex v in gProduct:
					if(getVertexId(v, gProduct) == 6)
						six = v;
				vm.insert(VertexMap::value_type(seven, six));
				// --- REPLACE setById(vm, 7, 6); STOP ---
			}
		} // end of Example TWO

		{ // Example THREE --- an actual enumeration algorithm
			// In essence, we iterate over all permutations of the vertices in the educt graph,
			// and for each of them we make a mapping to the vertices of the product graph.
			// Note, this does no 'alternating addition/removal' check for the edges, so you will very likely
			// get non-chemical reactions. This will in essence then be your main task:
			// make sure you will have only chemically valid rules as described on the webpage.
			// While it can be done by skipping permutations, a fast solution would probably not rely on simple
			// permutation generation.
			// For now we just disable the chemistry check as a hax to get something printed:
			std::cout << "WARNING: doChemistryCheck forced to false in " << __FILE__
			          << " at line " << (__LINE__ + 1) << "." << std::endl;
			doChemistryCheck = false; // this should be removed

			// Copy the educt verties into a vector so we can permute them.
			const auto vsEduct = vertices(gEduct);
			std::vector<UVertex> eductVertices(vsEduct.first, vsEduct.second);
			const auto productVertices = asRange(vertices(gEduct));

			// std::next_permutation uses sortedness to determine the last permuation, so start by sorting.
			std::sort(eductVertices.begin(), eductVertices.end());

			constexpr int limit = 10;
			// We will limit to 'limit' permutations, as otherwise the PDF might get too large.
			// Note that we also reject all rules where the atom type would change.
			int permutationCount = 0;
			do {
				VertexMap vm;
				bool valid = true;
				for(int i = 0; i != eductVertices.size(); i++) {
					// check if the molecule data for the atoms matches, otherwise we can't make a rule
					if(mod::graph::internal::getMolecule(eductVertices[i], lgEduct) !=
					   mod::graph::internal::getMolecule(productVertices[i], lgProduct)) {
						valid = false;
						break;
					}
					vm.insert(VertexMap::value_type(eductVertices[i], productVertices[i]));
				}
				if(!valid) continue;
				vertexMaps.push_back(std::move(vm));
				++permutationCount;
				if(permutationCount == limit) break;
			} while(std::next_permutation(eductVertices.begin(), eductVertices.end()));
		} // end of Example THREE
	} // END WORKING (roughly) HERE. 

	{ // debug stuff
		std::cout << "Maps:\n";
		for(const auto &vm : vertexMaps) {
			for(const auto &vt : vm) {
				std::cout << "\"" << mod::graph::internal::getMolecule(vt.left, lgEduct) << "\" "
				          << getVertexId(vt.left, gEduct)
				          << "\t<=>   "
				          << getVertexId(vt.right, gProduct) << " \""
				          << mod::graph::internal::getMolecule(vt.right, lgProduct) << "\"\n";
			}
			std::cout << '\n';
		}
	}

	std::vector<std::shared_ptr<mod::rule::Rule>> rules;
	for(const VertexMap &vm : vertexMaps) {
		auto r = createRule(vm, lgEduct, lgProduct, doChemistryCheck);
		// use r->setName(name); to give the rule a nicer name
		const auto iter = std::find_if(
				rules.begin(), rules.end(),
				[&r](const auto &rOther) {
					auto labelSettings = mod::LabelSettings(mod::LabelType::String,
					                                        mod::LabelRelation::Isomorphism);
					return 1 == r->isomorphism(rOther, 1, labelSettings);
				});
		if(iter != rules.end()) {
			std::cout << "Duplicate rule deleted\n";
		} else {
			rules.push_back(r);
		}
	}
	return rules;
}

#ifdef AS_PYTHON_EXTENSION
#include <boost/python.hpp>
#include <mod/Config.hpp>

namespace py = boost::python;

namespace {
// this can be used to make sure the extension and mod is using the same shared library

uintptr_t magicLibraryValue() {
	return (uintptr_t)&mod::getConfig();
}

} // namespace

BOOST_PYTHON_MODULE(pydoStuff) { // this macro argument is the name of the module, it must be the same as the .so file name.
	py::def("magicLibraryValue", &magicLibraryValue);

	// Change the string in the first argument to give the function another name in Python.
	// The second argument is the function pointer to the function above.
	py::def("doStuff", &doStuff);
}

#else // not AS_PYTHON_EXTENSION

int main(int argc, char **argv) {
	std::vector<std::shared_ptr<mod::graph::Graph> > educts, products;
	{ // do something else, e.g., take SMILES strings from argv
		std::shared_ptr<mod::graph::Graph> g1, g2;
		g1 = mod::graph::Graph::fromSMILES("OCC=O");
		g2 = mod::graph::Graph::fromSMILES("OC=CO");
		educts.push_back(g1);
		products.push_back(g2);
	}
	auto rules = doStuff(educts, products, true);
	for(auto r : rules) r->print();
	return 0;
}

#endif // AS_PYTHON_EXTENSION

//------------------------------------------------------------------------------
// Library stuff
//------------------------------------------------------------------------------

#include <mod/graph/Printer.hpp>
#include <mod/rule/Rule.hpp>
#include <mod/lib/Graph/Properties/String.hpp>
#include <mod/lib/Rules/Properties/Molecule.hpp>

std::shared_ptr<mod::rule::Rule> createRule(const VertexMap &vertexMap,
                                            const LUG &lgEduct, const LUG &lgProduct, bool doChemistryCheck) {
	using namespace mod;
	using namespace mod::lib;
	const auto &gEduct = mod::graph::internal::getGraph(lgEduct);
	const auto &gProduct = mod::graph::internal::getGraph(lgProduct);
	const auto n = num_vertices(gEduct);
	{ // error checking
		const auto n2 = num_vertices(gProduct);
		if(n != n2) {
			std::cout << "Different number of vertices: " << n << " and " << n2 << std::endl;
			assert(false);
			std::exit(1);
		}

		for(const auto &vt : vertexMap) {
			const auto vLeft = vt.left;
			const auto vRight = vt.right;
			const auto vLeftId = get(boost::vertex_index_t(), gEduct, vLeft);
			const auto vRightId = get(boost::vertex_index_t(), gProduct, vRight);
			if(vLeftId >= n) {
				std::cout << "Invalid left vertex, index is " << vLeftId << std::endl;
				assert(false);
				std::exit(1);
			}
			if(vRightId >= n) {
				std::cout << "Invalid right vertex, index is " << vRightId << std::endl;
				assert(false);
				std::exit(1);
			}
			const auto &lLabel = mod::graph::internal::getString(vLeft, lgEduct);
			const auto &rLabel = mod::graph::internal::getString(vRight, lgProduct);
			if(lLabel != rLabel) {
				std::cout << "Label mismatch: " << lLabel << " " << vLeftId << " <=> " << vRightId << " " << rLabel
				          << std::endl;
				assert(false);
				std::exit(1);
			}
		}
	} // end of error checking

	using CoreVertex = mod::lib::Rules::Vertex;
	using CoreEdge = mod::lib::Rules::Edge;
	using Membership = mod::lib::Rules::Membership;
	auto dpoRule = mod::rule::internal::makeLabelledRule();
	auto &core = mod::rule::internal::getGraph(dpoRule);
	dpoRule.pString = mod::rule::internal::makePropStringCore(core);
	auto &pStringCore = *dpoRule.pString;
	std::vector<CoreVertex> leftToCore(n, core.null_vertex()), rightToCore(n, core.null_vertex());

	// copy all matched vertices form the left
	for(const auto v : asRange(vertices(gEduct))) {
		const auto iter = vertexMap.left.find(v);
		if(iter == vertexMap.left.end()) continue;
		const auto vId = get(boost::vertex_index_t(), gEduct, v);
		const auto &label = mod::graph::internal::getString(v, lgEduct);
		const auto vCore = add_vertex(core);
		core[vCore].membership = Membership::Context;
		mod::rule::internal::add(pStringCore, vCore, label, label);
		leftToCore[vId] = vCore;
	}

	// link vertices from right side
	for(const auto v : asRange(vertices(gProduct))) {
		const auto iter = vertexMap.right.find(v);
		if(iter == vertexMap.right.end())continue;
		const auto vId = get(boost::vertex_index_t(), gProduct, v);
		const auto left = iter->second;
		const auto leftId = get(boost::vertex_index_t(), gEduct, left);
		const auto vCore = leftToCore[leftId];
		assert(core[vCore].membership == Membership::Context);
		rightToCore[vId] = vCore;
	}

	// copy left edges
	for(const auto e : asRange(edges(gEduct))) {
		const auto src = source(e, gEduct);
		const auto tar = target(e, gEduct);
		const auto srcCore = leftToCore[get(boost::vertex_index_t(), gEduct, src)];
		const auto tarCore = leftToCore[get(boost::vertex_index_t(), gEduct, tar)];
		if(srcCore == core.null_vertex()) continue;
		if(tarCore == core.null_vertex()) continue;
		const auto &label = mod::graph::internal::getString(e, lgEduct);
		const auto p = add_edge(srcCore, tarCore, core);
		assert(p.second);
		core[p.first].membership = Membership::Left;
		mod::rule::internal::add(pStringCore, p.first, label, "");
	}

	// copy right edges, or promote to context
	for(const auto e : asRange(edges(gProduct))) {
		const auto src = source(e, gProduct);
		const auto tar = target(e, gProduct);
		const auto srcCore = rightToCore[get(boost::vertex_index_t(), gProduct, src)];
		const auto tarCore = rightToCore[get(boost::vertex_index_t(), gProduct, tar)];
		if(srcCore == core.null_vertex()) continue;
		if(tarCore == core.null_vertex()) continue;
		auto p = edge(srcCore, tarCore, core);
		const auto &label = mod::graph::internal::getString(e, lgProduct);
		if(p.second) {
			core[p.first].membership = Membership::Context;
			mod::rule::internal::setRight(pStringCore, p.first, label);
		} else {
			p = add_edge(srcCore, tarCore, core);
			assert(p.second);
			core[p.first].membership = Membership::Right;
			mod::rule::internal::add(pStringCore, p.first, "", label);
		}
	}

	if(doChemistryCheck) {
		const auto error = [&]() {
			auto r = mod::rule::internal::makeRule(std::move(dpoRule));
			r->print();
			std::cout << "Run 'mod_post' to see rule." << std::endl;
			std::exit(1);
		};

		const auto molView = mod::rule::internal::makePropMoleculeCore(core, pStringCore);
		for(const auto v : asRange(vertices(core))) {
			unsigned int left = 0, right = 0;
			int bondChange = 0;
			for(const auto e : asRange(out_edges(v, core))) {
				const auto membership = core[e].membership;
				if(membership == Membership::Left) {
					left++;
					const auto bt = get(molView.getLeft(), e);
					bondChange -= bondValue(bt);
				} else if(membership == Membership::Right) {
					right++;
					const auto bt = get(molView.getRight(), e);
					bondChange += bondValue(bt);
				} else {
					const auto btLeft = get(molView.getLeft(), e);
					const auto btRight = get(molView.getRight(), e);
					if(btLeft == btRight) continue;
					left++;
					right++;
					bondChange -= bondValue(btLeft);
					bondChange += bondValue(btRight);
				}
			}
			if(bondChange != 0) {
				std::cout << "Non-zero bond change for vertex " << get(boost::vertex_index_t(), core, v) << ": "
				          << bondChange << "\n";
				std::cout << "Label: (" << get(pStringCore.getLeft(), v) << ", " << get(pStringCore.getRight(), v) << ")"
				          << std::endl;
				error();
			}
			if(left > 2) {
				std::cout << "Too many left side edges (" << left << ") for vertex "
				          << get(boost::vertex_index_t(), core, v) << ".\n";
				std::cout << "Label: (" << get(pStringCore.getLeft(), v) << ", " << get(pStringCore.getRight(), v) << ")"
				          << std::endl;
				error();
			}
			if(right > 2) {
				std::cout << "Too many right side edges (" << right << ") for vertex "
				          << get(boost::vertex_index_t(), core, v) << ".\n";
				std::cout << "Label: (" << get(pStringCore.getLeft(), v) << ", " << get(pStringCore.getRight(), v) << ")"
				          << std::endl;
				error();
			}
		}
	}

	return mod::rule::internal::makeRule(std::move(dpoRule));
}
