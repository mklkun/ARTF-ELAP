#include <sys/time.h>

#include <omp.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <deque>
#include <map>

#include <boost/multi_array.hpp>

#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/adjacency_iterator.hpp>

#include <boost/graph/make_connected.hpp>

#include <boost/graph/make_biconnected_planar.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/graph/random.hpp>

#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <boost/graph/grid_graph.hpp>
#include <boost/graph/copy.hpp>

#include <boost/graph/boyer_myrvold_planar_test.hpp>

#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/graph/small_world_generator.hpp>
#include <boost/graph/plod_generator.hpp>
#include <boost/random/linear_congruential.hpp>




using namespace std;
using namespace boost;



enum colors { WHITE, BLACK, RED ,GREY, GREEN };
enum sommets { A, B, C, D, E, F };



  typedef adjacency_matrix
    < undirectedS,
      property <vertex_index_t, int, property <vertex_degree_t, int > >,
      property <edge_index_t, int, property <edge_color_t, int > >
    >
    graphL;

  typedef adjacency_list
    < setS,
	  vecS,
	  undirectedS,
      property <vertex_index_t, int, property <vertex_degree_t, int > >,
      property <edge_index_t, int, property <edge_color_t, int > >
    >
    graphM;


    typedef graph_traits<graphM>::vertex_descriptor vertex_descriptorM;
    typedef graph_traits<graphM>::vertices_size_type vertices_size_typeM;
    typedef graph_traits<graphM>::vertex_iterator vertex_iteratorM;
    typedef graph_traits<graphM>::edge_descriptor edge_descriptorM;
    typedef graph_traits<graphM>::edges_size_type edges_size_typeM;
    typedef graph_traits<graphM>::edge_iterator edge_iteratorM;



    typedef boost::graph_traits<graphM> graph_traitsM;


namespace myNameSpace
{
    typedef boost::grid_graph < 2 > Grid;
    typedef boost::graph_traits < Grid > GridTraits;

    struct grid_to_graph_vertex_copier
    {
        typedef boost::property_map< Grid, boost::vertex_index_t>::type grid_vertex_index_map;
        typedef boost::property_map< ::graphM, boost::vertex_index_t>::type graph_vertex_index_map;

        const Grid& grid;
        grid_vertex_index_map grid_vertex_index;
        graph_vertex_index_map graph_vertex_index;

        grid_to_graph_vertex_copier(const Grid& grid_, graphM & graph) : grid(grid_), grid_vertex_index(get(boost::vertex_index_t(), grid_)), graph_vertex_index(get(boost::vertex_index_t(), graph))
        {
        }



    public:
        void operator() (GridTraits::vertex_descriptor grid_vertex, ::vertex_descriptorM graph_vertex) const
        {
            //size_t idx = get(grid_vertex_index, grid_vertex);
            //graph_vertex_index[graph_vertex] = idx;
        }
    };

    struct grid_to_graph_edge_copier
    {
        void operator() (GridTraits::edge_descriptor grid_edge, ::graph_traitsM::edge_descriptor graph_edge) const
        {
        }
    };
}




    bool extract_number(string & str, int & number)
    {
        string temp;

        unsigned int i=0;

        while ( i < str.size() )
        {
            if(isdigit(str[i]))
            {
                while (i < str.size() && isdigit(str[i]))
                {
                    temp += str[i];
                    i++;
                }
                istringstream stream(temp);
                stream >> number;

                if(i == str.size())
                    str.clear();
                else
                    str.erase( 0, i);

                return(true);
            }

            i++;
        }

        str.clear();
        return(false);
    }


    template <typename GRAPH>
    pair <edge_descriptorM, bool> get_edge(int index, GRAPH g)
    {
        edge_iteratorM ei, ei_end;
        pair <edge_descriptorM, bool> k;
        k.second = false;

        for (tuples::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
            if(get(edge_index, g, *ei) == index)
            {
                k.first = *ei;
                k.second = true;
            }
        return(k);
    }


    template <typename CYCLEIDS, typename CYCLES>
    bool appartient_a_non_corded(edge_descriptorM e, CYCLEIDS cycles_non_corded, CYCLES cycles)
    {
        for(unsigned i(0); i < cycles_non_corded.size(); i++)
        {
            unsigned temp = cycles_non_corded[i];
            for(unsigned j(0); j < cycles[temp].size(); j++)
                if(cycles[temp][j] == e)
                    return(true);
        }
        return(false);
    }


/////////////////////////////////////
//////       tri fusion        //////
/////////////////////////////////////


    template <typename V, typename GRAPH>
    int calculateDegreesSucc(V v, GRAPH g)
    {
        int deg(0);

        typename graph_traits< GRAPH >::out_edge_iterator out_i, out_end;

        for (boost::tuples::tie(out_i, out_end) = out_edges(v, g); out_i != out_end; ++out_i)
            deg += get(vertex_degree, g, target(*out_i, g));

        return(deg);
    }


    template <typename V, typename GRAPH>
    bool degresSuccSuperieurs(V v1, V v2, GRAPH g)
    {
        int deg1(0), deg2(0);

        deg1 = calculateDegreesSucc(v1, g);
        deg2 = calculateDegreesSucc(v2, g);


        if ( deg1 > deg2 )
            return(true);

        else if ( deg1 == deg2 )
        {
            typename graph_traits< GRAPH >::out_edge_iterator out_i, out_end;

            for (boost::tuples::tie(out_i, out_end) = out_edges(v1, g); out_i != out_end; ++out_i)
                deg1 += calculateDegreesSucc(target(*out_i, g), g);

            for (boost::tuples::tie(out_i, out_end) = out_edges(v2, g); out_i != out_end; ++out_i)
                deg2 += calculateDegreesSucc(target(*out_i, g), g);

            if ( deg1 > deg2 )
                return(true);
        }

        return(false);
    }


    template <typename V, typename GRAPH>
    void fusionerVertices(V & t,const int debut1,const int fin1,const int fin2, GRAPH g)
    {
        V t2;
        int debut2 = fin1+1;
        int compteur1 = debut1;
        int compteur2 = debut2;
        int i;

        //tableau2 = (int*)malloc((fin1-debut1+1)*sizeof(int));

        // copie des éléments du début de tableau
        for(i=debut1; i<=fin1; i++)
            t2.push_back(t[i]);

        // fusion des deux tableaux
        for(i=debut1; i<=fin2; i++)
        {
            if(compteur1==debut2) // éléments du 1er tableau tous utilisés
                break; // éléments tous classés


            if(compteur2==(fin2+1)) // éléments du 2nd tableau tous utilisés
            { // copie en fin de tableau des éléments du 1er sous tableau
                t[i] = t2[compteur1-debut1];
                compteur1++;
            }

            else if( get(vertex_degree, g, t2[compteur1-debut1]) > get(vertex_degree, g, t[compteur2])  )                                  //tableau2[compteur1-debut1]<tableau[compteur2])
            { // ajout d'1 élément du 1er sous tableau
                t[i] = t2[compteur1-debut1];
                compteur1++; // on avance ds le 1er sous tableau
            }

            else if( get(vertex_degree, g, t2[compteur1-debut1]) == get(vertex_degree, g, t[compteur2]) and degresSuccSuperieurs(t2[compteur1-debut1], t[compteur2], g) )                                  //tableau2[compteur1-debut1]<tableau[compteur2])
            {
                 // ajout d'1 élément du 1er sous tableau
                t[i] = t2[compteur1-debut1];
                compteur1++; // on avance ds le 1er sous tableau
            }

            else
            { // copie de l'élément à la suite du tableau
                t[i] = t[compteur2];
                compteur2++; // on avance ds le 2nd sous tableau
            }
        }
        //free(tableau2);
    }


    template <typename V, typename GRAPH>
    void triFusionAuxVertices(V & t, const int debut, const int fin, GRAPH g)
    {
        if(debut!=fin) // condition d'arrêt
        {
            int milieu = (debut+fin)/2;

            #pragma omp parallel sections
            {
                #pragma omp section
                triFusionAuxVertices(t, debut, milieu, g); // trie partie1
                #pragma omp section
                triFusionAuxVertices(t, milieu+1, fin, g); // trie partie2
            }
            fusionerVertices(t, debut, milieu, fin, g); // fusion des deux parties
        }
    }


    template <typename V, typename GRAPH>
    void triFusionVertices(V & t, const int longueur, GRAPH g)
    {
        if(longueur>0)
            triFusionAuxVertices(t, 0, longueur-1, g);
    }




    template <typename V, typename GRAPH, typename DEGREES>
    int calculateDegreesSuccR(V v, GRAPH g, DEGREES degrees)
    {
        int deg(0);

        typename graph_traits< GRAPH >::out_edge_iterator out_i, out_end;

        for (boost::tuples::tie(out_i, out_end) = out_edges(v, g); out_i != out_end; ++out_i)
            deg += degrees[target(*out_i, g)];

        return(deg);
    }


    template <typename V, typename GRAPH, typename DEGREES>
    bool degresSuccSuperieursR(V v1, V v2, GRAPH g, DEGREES degrees)
    {
        int deg1(0), deg2(0);

        deg1 = calculateDegreesSucc(v1, g);
        deg2 = calculateDegreesSucc(v2, g);


        if ( deg1 > deg2 )
            return(true);

        else if ( deg1 == deg2 )
        {
            typename graph_traits< GRAPH >::out_edge_iterator out_i, out_end;

            for (boost::tuples::tie(out_i, out_end) = out_edges(v1, g); out_i != out_end; ++out_i)
                deg1 += calculateDegreesSuccR(target(*out_i, g), g, degrees);

            for (boost::tuples::tie(out_i, out_end) = out_edges(v2, g); out_i != out_end; ++out_i)
                deg2 += calculateDegreesSuccR(target(*out_i, g), g, degrees);

            if ( deg1 > deg2 )
                return(true);
        }

        return(false);
    }


    template <typename V, typename GRAPH, typename DEGREES>
    void fusionerVerticesR(V & t,const int debut1,const int fin1,const int fin2, GRAPH g, DEGREES deg)
    {
        V t2;
        int debut2 = fin1+1;
        int compteur1 = debut1;
        int compteur2 = debut2;
        int i;

        //tableau2 = (int*)malloc((fin1-debut1+1)*sizeof(int));

        // copie des éléments du début de tableau
        for(i=debut1; i<=fin1; i++)
            t2.push_back(t[i]);


        // fusion des deux tableaux
        for(i=debut1; i<=fin2; i++)
        {
            if(compteur1==debut2) // éléments du 1er tableau tous utilisés
                break; // éléments tous classés

            else if(compteur2==(fin2+1)) // éléments du 2nd tableau tous utilisés
            { // copie en fin de tableau des éléments du 1er sous tableau
                t[i] = t2[compteur1-debut1];
                compteur1++;
            }

            else if( deg[t2[compteur1-debut1]] > deg[t[compteur2]]  )                                  //tableau2[compteur1-debut1]<tableau[compteur2])
            { // ajout d'1 élément du 1er sous tableau
                t[i] = t2[compteur1-debut1];
                compteur1++; // on avance ds le 1er sous tableau
            }

            else if( deg[t2[compteur1-debut1]] == deg[t[compteur2]] and degresSuccSuperieursR(t2[compteur1-debut1], t[compteur2], g, deg) )                                  //tableau2[compteur1-debut1]<tableau[compteur2])
            {
                 // ajout d'1 élément du 1er sous tableau
                t[i] = t2[compteur1-debut1];
                compteur1++; // on avance ds le 1er sous tableau
            }

            else
            { // copie de l'élément à la suite du tableau
                t[i] = t[compteur2];
                compteur2++; // on avance ds le 2nd sous tableau
            }
        }
        //free(tableau2);
    }


    template <typename V, typename GRAPH, typename DEGREES>
    void triFusionAuxVerticesR(V & t, const int debut, const int fin, GRAPH g, DEGREES deg)
    {
        if(debut!=fin) // condition d'arrêt
        {
            int milieu = (debut+fin)/2;

            #pragma omp parallel sections
            {
                #pragma omp section
                triFusionAuxVerticesR(t, debut, milieu, g, deg); // trie partie1
                #pragma omp section
                triFusionAuxVerticesR(t, milieu+1, fin, g, deg); // trie partie2
            }
            fusionerVerticesR(t, debut, milieu, fin, g, deg); // fusion des deux parties
        }
    }


    template <typename V, typename GRAPH, typename DEGREES>
    void triFusionVerticesR(V & t, const int longueur, GRAPH g, DEGREES deg)
    {
        if(longueur>0)
            triFusionAuxVerticesR(t, 0, longueur-1, g, deg);
    }




    template <typename I>
    void fusionerInts(I & t,const int debut1,const int fin1,const int fin2, const char c)
    {
        I t2;
        int debut2 = fin1+1;
        int compteur1 = debut1;
        int compteur2 = debut2;
        int i;

        //tableau2 = (int*)malloc((fin1-debut1+1)*sizeof(int));

        // copie des éléments du début de tableau
        for(i=debut1; i<=fin1; i++)
            t2.push_back(t[i]);


        // fusion des deux tableaux
        for(i=debut1; i<=fin2; i++)
        {
            if(compteur1==debut2) // éléments du 1er tableau tous utilisés
                break; // éléments tous classés

            else if(compteur2==(fin2+1)) // éléments du 2nd tableau tous utilisés
            { // copie en fin de tableau des éléments du 1er sous tableau
                t[i] = t2[compteur1-debut1];
                compteur1++;
            }

            else if( t2[compteur1-debut1] >= t[compteur2] && c == '>')                                  //tableau2[compteur1-debut1]<tableau[compteur2])
            { // ajout d'1 élément du 1er sous tableau
                t[i] = t2[compteur1-debut1];
                compteur1++; // on avance ds le 1er sous tableau
            }

            else if( t2[compteur1-debut1] <= t[compteur2] && c == '<')                                  //tableau2[compteur1-debut1]<tableau[compteur2])
            { // ajout d'1 élément du 1er sous tableau
                t[i] = t2[compteur1-debut1];
                compteur1++; // on avance ds le 1er sous tableau
            }

            else
            { // copie de l'élément à la suite du tableau
                t[i] = t[compteur2];
                compteur2++; // on avance ds le 2nd sous tableau
            }
        }
        //free(tableau2);
    }


    template <typename I>
    void triFusionAuxInts(I & t, const int debut, const int fin, const char c)
    {
        if(debut!=fin) // condition d'arrêt
        {
            int milieu = (debut+fin)/2;

            #pragma omp parallel sections
            {
                #pragma omp section
                triFusionAuxInts(t, debut, milieu, c); // trie partie1
                #pragma omp section
                triFusionAuxInts(t, milieu+1, fin, c); // trie partie2
            }
            fusionerInts(t, debut, milieu, fin, c); // fusion des deux parties
        }
    }


    template <typename I>
    void triFusionInts(I & t, const int longueur, const char c)
    {
        if(longueur>0)
            triFusionAuxInts(t, 0, longueur-1, c);
    }




    template <typename D>
    void fusionerDeques(D & t,const int debut1,const int fin1,const int fin2, const char c)
    {
        D t2;
        int debut2 = fin1+1;
        int compteur1 = debut1;
        int compteur2 = debut2;
        int i;

        //tableau2 = (int*)malloc((fin1-debut1+1)*sizeof(int));

        // copie des éléments du début de tableau
        for(i=debut1; i<=fin1; i++)
            t2.push_back(t[i]);


        // fusion des deux tableaux
        for(i=debut1; i<=fin2; i++)
        {
            if(compteur1==debut2) // éléments du 1er tableau tous utilisés
                break; // éléments tous classés

            else if(compteur2==(fin2+1)) // éléments du 2nd tableau tous utilisés
            { // copie en fin de tableau des éléments du 1er sous tableau
                t[i] = t2[compteur1-debut1];
                compteur1++;
            }

            else if( t2[compteur1-debut1].size() >= t[compteur2].size() && c == '>')                                  //tableau2[compteur1-debut1]<tableau[compteur2])
            { // ajout d'1 élément du 1er sous tableau
                t[i] = t2[compteur1-debut1];
                compteur1++; // on avance ds le 1er sous tableau
            }

            else if( t2[compteur1-debut1].size() <= t[compteur2].size() && c == '<')                                  //tableau2[compteur1-debut1]<tableau[compteur2])
            { // ajout d'1 élément du 1er sous tableau
                t[i] = t2[compteur1-debut1];
                compteur1++; // on avance ds le 1er sous tableau
            }

            else
            { // copie de l'élément à la suite du tableau
                t[i] = t[compteur2];
                compteur2++; // on avance ds le 2nd sous tableau
            }
        }
        //free(tableau2);
    }


    template <typename D>
    void triFusionAuxDeques(D & t, const int debut, const int fin, const char c)
    {
        if(debut!=fin) // condition d'arrêt
        {
            int milieu = (debut+fin)/2;

            #pragma omp parallel sections
            {
                #pragma omp section
                triFusionAuxDeques(t, debut, milieu, c); // trie partie1
                #pragma omp section
                triFusionAuxDeques(t, milieu+1, fin, c); // trie partie2
            }
            fusionerDeques(t, debut, milieu, fin, c); // fusion des deux parties
        }
    }


    template <typename D>
    void triFusionDeques(D & t, const int longueur, const char c)
    {
        if(longueur>0)
            triFusionAuxDeques(t, 0, longueur-1, c);
    }


////////////////////////////////////
////////////////////////////////////
////////////////////////////////////




    template <typename VD, typename EC, typename GRAPH, typename DED, typename DI, typename EI>
    void calculate_cycle(const VD s, const VD d, const EC colors, const GRAPH gm, DED & cycles, DI & e_cycles, const int c, const EI indexes)
    {
        typename graph_traits< GRAPH >::out_edge_iterator out_i, out_end;
        for (boost::tuples::tie(out_i, out_end) = out_edges(s, gm); out_i != out_end; ++out_i)
        {
            edge_descriptorM e(*out_i);

            if (colors[e] == GREEN)
            {
                if (depth_search_cycles(s, target(e, gm), d, colors, gm, cycles, e_cycles, c, indexes))
                {
                    cycles[c].push_back(e);
                    e_cycles[indexes[e]].push_back(c);
                }
            }
        }
    }


    template <typename VD, typename EC, typename GRAPH, typename DED, typename DI, typename EI>
    bool depth_search_cycles(const VD parent, const VD s, const VD d, const EC colors, const GRAPH gm, DED & cycles, DI & e_cycles, const int c, const EI indexes)
    {
        typename graph_traits< GRAPH >::out_edge_iterator out_i, out_end;
        if (s == d)
        {
            return(true);
        }
        for (boost::tuples::tie(out_i, out_end) = out_edges(s, gm); out_i != out_end; ++out_i)
        {
            edge_descriptorM e(*out_i);

            if (colors[e] == GREEN && target(e, gm) != parent)
            {
                if (depth_search_cycles(s, target(e, gm), d, colors, gm, cycles, e_cycles, c, indexes))
                {
                    cycles[c].push_back(e);
                    e_cycles[indexes[e]].push_back(c);
                    return(true);
                }
            }
        }
        return(false);
    }





    template < typename CYCLE, typename GRAPH >
    void combine_cycles(const CYCLE cycle_i, const CYCLE cycle_j, CYCLE & cycle_result, const GRAPH gm)
    {
        CYCLE temp_i(cycle_i), temp_j(cycle_j);
        deque < edge_descriptorM >::iterator pos;
        bool flag(false);

        for (unsigned i=0; i < cycle_i.size(); i++)
        {
            pos = find(temp_j.begin(), temp_j.end(), cycle_i[i]);
            if ( pos != temp_j.end())
            {
                temp_j.erase(pos);
                pos = find(temp_i.begin(), temp_i.end(), cycle_i[i]);
                temp_i.erase(pos);
                flag = true;
            }
        }
        if(flag)
        {
            temp_j.insert(temp_j.begin(), temp_i.begin(), temp_i.end());
            cycle_result.insert(cycle_result.begin(), temp_j.begin(), temp_j.end());

            edge_descriptorM e(temp_j.front());
            temp_j.pop_front();

            if(!is_cycle(source(e, gm), target(e, gm), temp_j, gm))
                cycle_result.clear();
        }
    }


    template < typename VERTEX, typename CYCLE , typename GRAPH>
    bool is_cycle(VERTEX s, VERTEX d, CYCLE & cycle, const GRAPH gm)
    {
        for (unsigned i=0; i < cycle.size(); i++)
        {
            if(source(cycle[i], gm) == s)
            {
                if(target(cycle[i], gm) == d)
                    return(true);
                VERTEX temp(target(cycle[i], gm));
                cycle.erase( cycle.begin() + i );
                return(is_cycle(temp, d, cycle, gm));
            }
            else if(source(cycle[i], gm) == d)
            {
                if(target(cycle[i], gm) == s)
                    return(true);
                VERTEX temp(target(cycle[i], gm));
                cycle.erase( cycle.begin() + i );
                return(is_cycle(s, temp, cycle, gm));
            }
            else if(target(cycle[i], gm) == s)
            {
                VERTEX temp(source(cycle[i], gm));
                cycle.erase( cycle.begin() + i );
                return(is_cycle(temp, d, cycle, gm));
            }
            else if(target(cycle[i], gm) == d)
            {
                VERTEX temp(source(cycle[i], gm));
                cycle.erase( cycle.begin() + i );
                return(is_cycle(s, temp, cycle, gm));
            }
        }
        return (false);
    }





int main()
{

    // ProcessTime example
    struct timeval startTime;
    struct timeval endTime;
    // get the current time
    // - NULL because we don't care about time zone
    gettimeofday(&startTime, NULL);

    double tS, tE;

    graphM gm;


    ofstream fichier_res("resultats.txt", ios::out | ios::trunc);


// Remplissage du Graphe
{

/*    boost::array < size_t, 2 > lengths = { { 20, 20 } };
    myNameSpace::Grid grid(lengths, false);


    boost::copy_graph(grid, gm, boost::vertex_copy(myNameSpace::grid_to_graph_vertex_copier(grid, gm)).edge_copy(myNameSpace::grid_to_graph_edge_copier()));
    cout << endl;  */


/*	vertices_size_typeM X(30);
    vertices_size_typeM Y(250);

    mt19937 rng(time(0));
    generate_random_graph( gm, X, Y, rng, false, false);


    std::ofstream filename("graphes/Random/random_test.dot");
    write_graphviz(filename,gm);  */


/*    vertices_size_typeM X(20);

    typedef boost::erdos_renyi_iterator<boost::minstd_rand, graphM> ERGen;
    minstd_rand gen;
    // Create graph with 100 nodes and edges with probability 0.05
    graphM gm(ERGen(gen, X, 0.7), ERGen(), X);   */




/*    vertices_size_typeM X(40);

    typedef plod_iterator<minstd_rand, graphM> SFGen;
    minstd_rand gen(time(0));
    // Create graph with 100 nodes
    graphM gm(SFGen(gen, X, 2.2, 5000), SFGen(), X);

    //std::ofstream filename("graphes/Random/random_10_02.dot");
    //write_graphviz(filename,gm);  */




   /* add_edge(A, B, gm);
    add_edge(A, C, gm);
    add_edge(D, B, gm);
    add_edge(D, E, gm);
    add_edge(F, C, gm);
    add_edge(F, E, gm);
    add_edge(E, B, gm);
    add_edge(C, B, gm);
    add_edge(C, E, gm);  exemple traingle benchmark avec 6 noeuds */



    //fstream fichier("test fichier.txt");



    //ifstream fichier("graphes/Old Random/random_10_20.dot", ios::in);
    //ifstream fichier("graphes/Graph_test/tunis.dot", ios::in);  // on ouvre en lecture
    ifstream fichier("graphes/Grids/grid_15_15.dot", ios::in);
    //ifstream fichier("graphes/Random/random_40_02.dot", ios::in);


    int number;
    vector <int> temp;

    vector < vector <int> > v;

    int num(0);

    if(fichier)  // si l'ouverture a fonctionné
    {
        string ligne;  // déclaration d'une chaîne qui contiendra la ligne lue
        while(getline(fichier, ligne))  // tant que l'on peut mettre la ligne dans "contenu"
        {
            while ( ! ligne.empty())
                if (extract_number(ligne, number))
                    temp.push_back(number);

            if (temp.size() == 1)
                num++;
            else if (temp.size() == 2 && temp[0]!=temp[1])
                v.push_back(vector <int> (temp));

            temp.clear();
        }

        fichier.close();
    }
    else
        return(10);


    for (unsigned i=0; i < v.size(); i++)
        add_edge(v[i][0], v[i][1], gm);

}


//graphe pre-traitments
{
    //make_connected(gm);

    //Test for planarity; compute the planar embedding as a side-effect
/*    vector< vector< edge_descriptorM > > embedding(num_vertices(gm));

    if (boyer_myrvold_planarity_test(boyer_myrvold_params::graph = gm, boyer_myrvold_params::embedding = &embedding[0]))
        make_biconnected_planar(gm, &embedding[0]);
    else
        cout << "Input graph is not planar" << endl;*/
}


    write_graphviz(std::cout, gm);

    //std::ofstream filename("graphes/Random/random_test.dot");
    //write_graphviz(filename,gm);



//vertices's initializations
{

    vertex_iteratorM ui, ui_end;

    //Creating a property_map with the degrees of the degrees of each vertex
    property_map<graphM,vertex_degree_t>::type deg = get(vertex_degree, gm);

    for (boost::tuples::tie(ui, ui_end) = vertices(gm); ui != ui_end; ++ui)
        deg[*ui] = out_degree(*ui, gm);


/*
    //Creating a property_map with our external properties of each vertex
    property_map<graphM,vertex_all_t>::type allV = get(vertex_all, gm);
    int k;
    for (k =0, boost::tuples::tie(ui, ui_end) = vertices(gm); ui != ui_end; ++ui, k++)
    {
		allV[*ui].v_flag = false;
		allV[*ui].v_index = k;
    }
*/

}


    property_map<graphM, edge_index_t>::type e_index = get(edge_index, gm);
    property_map<graphM, edge_color_t>::type e_color = get(edge_color, gm);


//edges's initializations
{
    edge_iteratorM ui, ui_end;

/*    property_map<graphM,edge_all_t>::type allE = get(edge_all, gm);
    int k;
    for (k =0, boost::tuples::tie(ui, ui_end) = edges(gm); ui != ui_end; ++ui, k++)
    {
		allE[*ui].e_flag = false;
		allE[*ui].e_index = k;
    }
*/

    edges_size_typeM edge_count = 0;

    for(boost::tuples::tie(ui, ui_end) = edges(gm); ui != ui_end; ++ui, edge_count++)
    {
        e_index[*ui] = edge_count;
        e_color[*ui] = BLACK;
    }

}


    deque < vertex_descriptorM > g_vertices;


//getting the vertices
{


    for (vertices_size_typeM s(0) ; s < num_vertices(gm) ; s++)
        g_vertices.push_back(vertex(s,gm));

    /*
    for (boost::tuples::tie(ui, ui_end) = vertices(gm); ui != ui_end; ++ui)
        g_vertices.push_back(*ui);
    //g_vertices.push_back(*ui_end);
    */

}


    triFusionVertices(g_vertices, g_vertices.size(), gm);

//affichage des sommets tries
{
/*    vertices_size_typeM s(0);

    for (s=0 ; s < g_vertices.size() ; s++)
        cout << "sommet " << s << " est d index : " << get(vertex_index, gm, g_vertices[s])  << " est de degre : " << get(vertex_degree, gm, g_vertices[s])  << endl;
    cout << endl;
*/
}

    cout << " le graphe contient " << num_vertices(gm) << " noeuds et " << num_edges(gm) << " aretes." << endl << endl;

    fichier_res << " le graphe contient " << num_vertices(gm) << " noeuds et " << num_edges(gm) << " aretes." << endl << endl;



    deque < edge_descriptorM > spanning;
    deque < edge_descriptorM > chords;

    //restart time
    gettimeofday(&startTime, NULL);


//calcul de l arbre de recouvrement (la solution initiale)
{
    deque < vertex_descriptorM > spanned_vertices;
	deque < unsigned > spanning_ids;
    deque < vertex_descriptorM > g_vertices_temp(g_vertices);


    graph_traits<graphM>::out_edge_iterator out_i, out_end;

    while ( g_vertices_temp.size() > 0 )
    {
        for (boost::tuples::tie(out_i, out_end) = out_edges(g_vertices_temp.front(), gm); out_i != out_end; ++out_i)
        {
            deque < vertex_descriptorM >::iterator pos1 = find(spanned_vertices.begin(), spanned_vertices.end(), g_vertices_temp.front());// or source(*out_i,gm));
            deque < vertex_descriptorM >::iterator pos2 = find(spanned_vertices.begin(), spanned_vertices.end(), target(*out_i,gm));

            if ( pos2 == spanned_vertices.end() )   //pos1 == spanned_vertices.end() ||
            {
				if ( pos1 == spanned_vertices.end() )
				{
					spanning.push_back(*out_i);
					e_color[*out_i] = GREEN;

                    spanned_vertices.push_back(target(*out_i,gm));
					spanning_ids.push_back(*max_element(spanning_ids.begin(),spanning_ids.end())+1);

                    spanned_vertices.push_back(g_vertices_temp.front());
					spanning_ids.push_back(spanning_ids.back());
				}
				else
				{
					spanning.push_back(*out_i);
					e_color[*out_i] = GREEN;

                    spanned_vertices.push_back(target(*out_i,gm));
					spanning_ids.push_back(spanning_ids[pos1 - spanned_vertices.begin()]);
				}
            }
            else if ( pos1 == spanned_vertices.end() )
            {
                spanning.push_back(*out_i);
                e_color[*out_i] = GREEN;

                spanned_vertices.push_back(g_vertices_temp.front());
				spanning_ids.push_back(spanning_ids[pos2 - spanned_vertices.begin()]);
            }
			else
			{
				if ( spanning_ids[pos1 - spanned_vertices.begin()] < spanning_ids[pos2 - spanned_vertices.begin()] )
				{
					spanning.push_back(*out_i);
					e_color[*out_i] = GREEN;
					for (unsigned c=0, temp=spanning_ids[pos2 - spanned_vertices.begin()]; c < spanning_ids.size(); c++)
						if (spanning_ids[c] == temp)
							spanning_ids[c] = spanning_ids[pos1 - spanned_vertices.begin()];
				}
				else if ( spanning_ids[pos1 - spanned_vertices.begin()] > spanning_ids[pos2 - spanned_vertices.begin()] )
				{
					spanning.push_back(*out_i);
					e_color[*out_i] = GREEN;
					for (unsigned c=0, temp=spanning_ids[pos1 - spanned_vertices.begin()]; c < spanning_ids.size(); c++)
						if (spanning_ids[c] == temp)
							spanning_ids[c] = spanning_ids[pos2 - spanned_vertices.begin()];
				}
				else
				{
				    if( e_color[*out_i] == BLACK)
				    {
				        chords.push_back(*out_i);
				        e_color[*out_i] = GREY;
				    }
				}
			}
        }
        g_vertices_temp.pop_front();
    }
}


//affichage des aretes formant l'arbre de recouvrement ainsi que les cordes
{
/*    for (unsigned i =0 ; i < spanning.size() ; i++)
        cout << "Branche " << i << " est d index : " << get(edge_index, gm, spanning[i])  << " entre les noeuds " << get(vertex_index, gm, source(spanning[i], gm)) << " et " << get(vertex_index, gm, target(spanning[i], gm)) <<endl;
    cout << endl;

    for (unsigned i =0 ; i < chords.size() ; i++)
        cout << "Corde " << i << " est d index : " << get(edge_index, gm, chords[i])  << " entre les noeuds " << get(vertex_index, gm, source(chords[i], gm)) << " et " << get(vertex_index, gm, target(chords[i], gm)) <<endl;
    cout << endl;  */
}


    deque < deque < edge_descriptorM > > cycles;
    deque < deque < int > > e_cycles;

    int tailleBCF(0);


//traitement d'initialisation des structures cycles et e_cycles
{
    for(unsigned i=0; i<chords.size(); ++i)
        cycles.push_back( deque < edge_descriptorM > () );

    for(unsigned i=0; i<num_edges(gm); ++i)
        e_cycles.push_back( deque <int> () );
}


//calcul des cycles
{
    #pragma omp parallel for
    for (unsigned c=0; c < chords.size(); c++)
    {
        cycles[c].push_back(chords[c]);
        e_cycles[get(edge_index, gm, chords[c])].push_back(c);
        calculate_cycle(source(chords[c], gm), target(chords[c], gm), e_color, gm, cycles, e_cycles, c, e_index);

        tailleBCF += cycles[c].size();
    }

    cout << "La taille de la BCF est de " << tailleBCF << endl << endl;

    fichier_res << tailleBCF << endl;
}


    gettimeofday(&endTime, NULL);
    // calculate time in microseconds
    tS = startTime.tv_sec*1000000 + (startTime.tv_usec);
    tE = endTime.tv_sec*1000000  + (endTime.tv_usec);
    std::cout << "Le calcul de solution initiale a dure : " << (tE - tS)/ 1000000.0L << endl << endl;

    fichier_res << (tE - tS)/ 1000000.0L << endl;

//fixing the gap between cycles and indexes
{
    for (unsigned i=0; i < cycles.size(); i++)
        for (unsigned j=0; j < cycles[i].size(); j++)
            cycles[i][j] = edge(source(cycles[i][j], gm), target(cycles[i][j], gm), gm).first;
}


//verifier si la BC est fondamentale
{

    deque < deque < int > > e_cordes;

    //extraction des cordes probables
    for (unsigned i=0; i < e_cycles.size(); i++)
        if (e_cycles[i].size() == 1)
        {
            e_cordes.push_back(deque < int > ());
            e_cordes.back().push_back(i);
            e_cordes.back().push_back(e_cycles[i].front());
        }


    for (unsigned i=0; i < e_cordes.size(); i++)
    {
        unsigned j(i+1);
        while (j < e_cordes.size())
        {
            if (e_cordes[i][1] == e_cordes[j][1])
                e_cordes.erase(e_cordes.begin() + j);
            else
                j++;
        }
    }


    if (e_cordes.size() < cycles.size())
    {
        cout << " La BC n'est pas fondamentale !!! " << endl;
        cout << "Car " << cycles.size() - e_cordes.size() << " cycles ne possedent pas de cordes !" << endl;
    }
    else
        cout << "Le BC est FONDAMENTALE !!! " << endl;

    //for (unsigned i=0; i < e_cordes.size(); i++)
        //cout << "l'arete " << e_cordes[i][0] << " appartient aux cycles : " << e_cordes[i][1] << endl;
}


//Calcul de la matrice d'adjacence des cycles
{

    // Create a 2D array
    typedef boost::multi_array<int, 2> matrix;
    //typedef matrix::index index;   //pour le parcours

    matrix mac(boost::extents[cycles.size()][cycles.size()]);

    for (unsigned i=0; i < e_cycles.size(); i++)
        for (unsigned j=0; j < e_cycles[i].size(); j++)
            for (unsigned k=j+1; k < e_cycles[i].size(); k++)
                mac[e_cycles[i][j]][e_cycles[i][k]] = 1;

    //Affichage de la MAC
    {
     /*   cout << endl;
        for (unsigned i=0; i < mac.size(); i++)
        {
            for (unsigned j=0; j < mac[i].size(); j++)
                cout << mac[i][j] << " ";
            cout << endl;
        } */
    }

    int numberOf1(0);

    //Calcul du nombre d'elements non nuls
    {
        for (unsigned i=0; i < mac.size(); i++)
            for (unsigned j=0; j < mac[i].size(); j++)
                if (mac[i][j]==1)
                    numberOf1++;
    }

    cout << endl << "Le nombre d'elements non nul est egal a " << numberOf1 << endl << endl;
    fichier_res << numberOf1 << endl << endl << endl << endl << endl;

    //verifier si triangulaire superieure
    {
        bool flag(true);
        for (unsigned i=0; i < mac.size(); i++)
            for (unsigned j=0; j < i; j++)
                if (mac[i][j]==1)
                    flag=false;
        if(flag)
            cout << "la matrice est triangulaire superieure." << endl << endl;
        else
            cout << "la matrice n'est pas triangulaire superieure." << endl << endl;
    }

}


//trier les cycles (cycles ou e_cyclces)
{
     //trier les cycles
     triFusionDeques(cycles, cycles.size(), '<');

     //trier les aretes selon leurs degres relatifs (a leurs appartenance aux cycles)
 //    triFusionDeques(e_cycles, e_cycles.size(), '<');     ne jamais le faire car l'indice i de e_cycle est la 'i'eme l'arete
}


//affichage des cycles et en fonction de ces derniers, et en fonction des aretes
{
/*   for (unsigned i=0; i < cycles.size(); i++)
    {
        cout << "le cycle " << i << " est de taille " << cycles[i].size() << " et contient les aretes :" << endl;
        for (unsigned j=0; j < cycles[i].size(); j++)
        {
            cout << "    l'arete " << get(edge_index, gm, cycles[i][j]) << " reliant les sommets " << get(vertex_index, gm, source(cycles[i][j], gm)) << " et " << get(vertex_index, gm, target(cycles[i][j], gm)) << endl;
        }
    }
    cout << endl;


    for (unsigned i=0; i < e_cycles.size(); i++)
    {
        cout << "l'arete " << i << " appartient aux cycles :" << endl;
        for (unsigned j=0; j < e_cycles[i].size(); j++)
        {
            cout << "    le cycle " << e_cycles[i][j] << endl;
        }
    }
    cout << endl; */

}


    //restart time
    gettimeofday(&startTime, NULL);

    int tailleBC(tailleBCF);


//heuristique iterative d'amelioration
{
    deque < edge_descriptorM > cycle_temp;
    int taillePrec(0),nbit(0);

    while( taillePrec-tailleBC )
    {
        nbit++;
        taillePrec = tailleBC;
        for (unsigned i=0; i < cycles.size(); i++)
        {
            #pragma omp parallel for private(cycle_temp)
            for (unsigned k=cycles.size(); k > i; k--)
            {
                unsigned j(k-1);
                //cout << " le nombre des cycles est " << cycles.size() << endl;
                //cout << "  i est egale a " << i << endl;
                //cout << "  j est egale a " << j << endl;
                if(!cycle_temp.empty())
                    cycle_temp.clear();
                if(i != j)
                {
                    combine_cycles(cycles[i], cycles[j], cycle_temp, gm);
                    if (cycle_temp.size()>0)
                    {
                        if ( cycles[i].size() < cycles[j].size() )
                        {
                            if ( cycle_temp.size() < cycles[i].size() )
                            {
                                // avec i ou avec j
                                tailleBC = ( tailleBC - cycles[j].size() ) + cycle_temp.size();
                                cycles[j].clear();
                                cycles[j].insert(cycles[j].begin(), cycle_temp.begin(), cycle_temp.end());
                            }
                            else if ( cycle_temp.size() < cycles[j].size() )
                            {
                                // changement avec j
                                tailleBC = ( tailleBC - cycles[j].size() ) + cycle_temp.size();
                                cycles[j].clear();
                                cycles[j].insert(cycles[j].begin(), cycle_temp.begin(), cycle_temp.end());
                            }
                        }
                        else
                        {
                            if ( cycle_temp.size() < cycles[j].size() )
                            {
                                // avec j ou avec i
                                tailleBC = ( tailleBC - cycles[i].size() ) + cycle_temp.size();
                                cycles[i].clear();
                                cycles[i].insert(cycles[i].begin(), cycle_temp.begin(), cycle_temp.end());
                            }
                            else if ( cycle_temp.size() < cycles[i].size() )
                            {
                                // changement avec i
                                tailleBC = ( tailleBC - cycles[i].size() ) + cycle_temp.size();
                                cycles[i].clear();
                                cycles[i].insert(cycles[i].begin(), cycle_temp.begin(), cycle_temp.end());
                            }
                        }
                    }
                }
            }
        }
        cout << " iteration : " << nbit << " taille precedente de la BC " << taillePrec << " nouvelle taille BC " << tailleBC<< endl;
    }



    cout << "La taille de la BCF initiale est de " << tailleBCF << endl << endl;

    cout << "La taille de la BC (F ou NF) apres l'HIA est de " << tailleBC << endl << endl;

    cout << "Il a fallut " << nbit << " iterations." << endl << endl;


    fichier_res << tailleBC << endl;
    fichier_res << nbit << endl;

}


    // get the end time
    gettimeofday(&endTime, NULL);
    // calculate time in microseconds
    tS = startTime.tv_sec*1000000 + (startTime.tv_usec);
    tE = endTime.tv_sec*1000000  + (endTime.tv_usec);
    std::cout << "L'heuristique iterative d'amelioration a dure : " << (tE - tS) / 1000000.0L << endl << endl;

    fichier_res << (tE - tS) / 1000000.0L << endl;


//recalculer e_cycles en fonction de la nouvelle base de cycles
{
    for (unsigned i=0; i < e_cycles.size(); i++)
        e_cycles[i].clear();

    for (unsigned i=0; i < cycles.size(); i++)
        for (unsigned j=0; j < cycles[i].size(); j++)
            e_cycles[get(edge_index, gm, cycles[i][j])].push_back(i);
}


//affichage des cycles et en fonction de ces derniers, et en fonction des aretes
{
/*   for (unsigned i=0; i < cycles.size(); i++)
    {
        cout << "le cycle " << i << " est de taille " << cycles[i].size() << " et contient les aretes :" << endl;
        for (unsigned j=0; j < cycles[i].size(); j++)
        {
            cout << "    l'arete " << get(edge_index, gm, cycles[i][j]) << " reliant les sommets " << get(vertex_index, gm, source(cycles[i][j], gm)) << " et " << get(vertex_index, gm, target(cycles[i][j], gm)) << endl;
        }
    }
    cout << endl;


    for (unsigned i=0; i < e_cycles.size(); i++)
    {
        cout << "l'arete " << i << " appartient aux cycles :" << endl;
        for (unsigned j=0; j < e_cycles[i].size(); j++)
        {
            cout << "    le cycle " << e_cycles[i][j] << endl;
        }
    }
    cout << endl;  */

}


    deque < deque < int > > e_cordes;


//verifier si la BC est fondamentale
{

    //edge_descriptorM e(get_edge(index, gm).first);

    //extraction des cordes probables
    for (unsigned i=0; i < e_cycles.size(); i++)
        if (e_cycles[i].size() == 1)
        {
            e_cordes.push_back(deque < int > ());
            e_cordes.back().push_back(i);
            e_cordes.back().push_back(e_cycles[i].front());
        }


    for (unsigned i=0; i < e_cordes.size(); i++)
    {
        unsigned j(i+1);
        while (j < e_cordes.size())
        {
            if (e_cordes[i][1] == e_cordes[j][1])
                e_cordes.erase(e_cordes.begin() + j);
            else
                j++;
        }
    }


    if (e_cordes.size() < cycles.size())
    {
        cout << " La BC n'est pas fondamentale !!! " << endl;
        cout << "Car " << cycles.size() - e_cordes.size() << " cycles ne possedent pas de cordes !" << endl;

        fichier_res << "NF(" << cycles.size() - e_cordes.size() << ")" << endl;
    }
    else
    {
        cout << "Le BC est FONDAMENTALE !!! " << endl;

        fichier_res << "F" << endl;
    }
    //for (unsigned i=0; i < e_cordes.size(); i++)
        //cout << "l'arete " << e_cordes[i][0] << " appartient aux cycles : " << e_cordes[i][1] << endl;
}


//Calcul de la matrice d'adjacence des cycles
{

    // Create a 2D array
    typedef boost::multi_array<int, 2> matrix;
    //typedef matrix::index index;   //pour le parcours

    matrix mac(boost::extents[cycles.size()][cycles.size()]);

    for (unsigned i=0; i < e_cycles.size(); i++)
        for (unsigned j=0; j < e_cycles[i].size(); j++)
            for (unsigned k=j+1; k < e_cycles[i].size(); k++)
                mac[e_cycles[i][j]][e_cycles[i][k]] = 1;

    //Affichage de la MAC
    {
     /*   cout << endl;
        for (unsigned i=0; i < mac.size(); i++)
        {
            for (unsigned j=0; j < mac[i].size(); j++)
                cout << mac[i][j] << " ";
            cout << endl;
        } */
    }

    int numberOf1(0);

    //Calcul du nombre d'elements non nuls
    {
        for (unsigned i=0; i < mac.size(); i++)
            for (unsigned j=0; j < mac[i].size(); j++)
                if (mac[i][j]==1)
                    numberOf1++;
    }

    cout << endl << "Le nombre d'elements non nul est egal a " << numberOf1 << endl << endl;

    fichier_res << numberOf1 << endl << endl << endl << endl << endl;


    //verifier si triangulaire superieure
    {
        bool flag(true);
        for (unsigned i=0; i < mac.size(); i++)
            for (unsigned j=0; j < i; j++)
                if (mac[i][j]==1)
                    flag=false;
        if(flag)
            cout << "la matrice est triangulaire superieure." << endl << endl;
        else
            cout << "la matrice n'est pas triangulaire superieure." << endl << endl;
    }

}


//verifier si le graphe est planaire
{
    //Test for planarity; compute the planar embedding as a side-effect
    typedef vector < edge_descriptorM > vec_e;
    vector < vec_e > embedding( num_vertices(gm) );
    if ( boyer_myrvold_planarity_test ( boyer_myrvold_params::graph = gm, boyer_myrvold_params::embedding = &embedding[0] ) )
    {
        cout << "Input graph is planar" << endl;

        fichier_res << "P"  << endl << endl << endl << endl;
    }
    else
    {
        cout << "Input graph is not planar" << endl;

        fichier_res << "NP"  << endl << endl << endl << endl;
    }
    cout << endl;
}


    //restart time
    gettimeofday(&startTime, NULL);


//rendre la base fondamentale
{
    deque < vertex_descriptorM > spanned_vertices;
    deque < unsigned > spanning_ids;
    deque < vertex_descriptorM > g_vertices_temp;

    //coloring all edges in RED
    {
        edge_iteratorM ei, ei_end;
        for(boost::tuples::tie(ei, ei_end) = edges(gm); ei != ei_end; ++ei)
            e_color[*ei] = RED;
    }

    //affichage des couleurs
    {
        /*edge_iteratorM ei, ei_end;
        for(boost::tuples::tie(ei, ei_end) = edges(gm); ei != ei_end; ++ei)
        {
            cout << "arete d index : " << get(edge_index, gm, *ei)  << " entre les noeuds " << get(vertex_index, gm, source(*ei, gm)) << " et " << get(vertex_index, gm, target(*ei, gm));
            if(e_color[*ei] == BLACK)
                cout << " est de couleur BLACK" << endl;
            if(e_color[*ei] == GREY)
                cout << " est de couleur GREY" << endl;
            if(e_color[*ei] == GREEN)
                cout << " est de couleur GREEN" << endl;
            if(e_color[*ei] == RED)
                cout << " est de couleur RED" << endl;
        }
        cout << endl;  */
    }

    deque < unsigned > cycles_corded;
    deque < unsigned > cycles_non_corded;

    chords.clear();

    //calculating cycles_corded and cycles_non_corded
    {
        for(unsigned i=0; i < e_cordes.size(); i++)
        {
            edge_descriptorM chord(get_edge(e_cordes[i][0], gm).first);
            chord = edge(source(chord,gm),target(chord,gm),gm).first;
            e_color[chord] = GREY;
            chords.push_back(chord);
            cycles_corded.push_back(e_cordes[i][1]);
        }


        triFusionInts(cycles_corded,cycles_corded.size(), '<');

        unsigned temp(cycles_corded.size()-1);

        for(unsigned i=0; i < cycles_corded[0]; i++)
            cycles_non_corded.push_back(i);

        for (unsigned i=0; i < temp; i++)
            for (unsigned j(cycles_corded[i]+1);j < cycles_corded[i+1];j++)
                cycles_non_corded.push_back(j);

        for(unsigned i=(cycles_corded.back()+1); i < cycles.size(); i++)
            cycles_non_corded.push_back(i);
    }

    //affichage des corded and non_corded
    {
       /* for(unsigned i=0; i < cycles_corded.size(); i++)
            cout << " corded cycle num " << cycles_corded[i] << endl;
        cout << " number of corded " << cycles_corded.size() << endl << endl;

        for(unsigned i=0; i < cycles_non_corded.size(); i++)
            cout << " non corded cycle num " << cycles_non_corded[i] << endl;
        cout << " number of non corded " << cycles_non_corded.size() << endl << endl; */
    }

    //adjusting non corded cycles
    {
        for (unsigned i=0; i < cycles_non_corded.size(); i++)
        {
            for (unsigned j=0; j < cycles[cycles_non_corded[i]].size(); j++)
            {
                edge_descriptorM e(cycles[cycles_non_corded[i]][j]);
                e_color[e] = BLACK;

                deque < vertex_descriptorM >::iterator pos ;

                pos = find(g_vertices_temp.begin(), g_vertices_temp.end(), source(e,gm));
                if ( pos == g_vertices_temp.end() )
                    g_vertices_temp.push_back(source(e,gm));

                pos = find(g_vertices_temp.begin(), g_vertices_temp.end(), target(e,gm));
                if ( pos == g_vertices_temp.end() )
                    g_vertices_temp.push_back(target(e,gm));
            }
        }
    }

    //affichage des couleurs
    {
      /*  edge_iteratorM ei, ei_end;
        for(boost::tuples::tie(ei, ei_end) = edges(gm); ei != ei_end; ++ei)
        {
            cout << "arete d index : " << get(edge_index, gm, *ei)  << " entre les noeuds " << get(vertex_index, gm, source(*ei, gm)) << " et " << get(vertex_index, gm, target(*ei, gm));
            if(e_color[*ei] == BLACK)
                cout << " est de couleur BLACK" << endl;
            if(e_color[*ei] == GREY)
                cout << " est de couleur GREY" << endl;
            if(e_color[*ei] == GREEN)
                cout << " est de couleur GREEN" << endl;
            if(e_color[*ei] == RED)
                cout << " est de couleur RED" << endl;
        }
        cout << endl; */
    }

    spanning.clear();

    //adjusting corded cycles
    {
        for(unsigned i=0; i < cycles_corded.size(); i++)
        {
            for(unsigned j=0; j < cycles[cycles_corded[i]].size(); j++)
            {
                edge_descriptorM e(cycles[cycles_corded[i]][j]);

                e = edge(source(e,gm),target(e,gm),gm).first;

                if ( e_color[e] != GREY and !appartient_a_non_corded(e, cycles_non_corded, cycles))
                {

                    deque < vertex_descriptorM >::iterator pos1 = find(spanned_vertices.begin(), spanned_vertices.end(), source(e,gm));
                    deque < vertex_descriptorM >::iterator pos2 = find(spanned_vertices.begin(), spanned_vertices.end(), target(e,gm));

                    if ( pos2 == spanned_vertices.end() )   //pos1 == spanned_vertices.end() ||
                    {
                        if ( pos1 == spanned_vertices.end() )
                        {
                            spanning.push_back(e);
					        e_color[e] = GREEN;

                            spanned_vertices.push_back(target(e,gm));
                            if(max_element(spanning_ids.begin(),spanning_ids.end()) == spanning_ids.end())
                                spanning_ids.push_back(0);
                            else
                                spanning_ids.push_back(*max_element(spanning_ids.begin(),spanning_ids.end())+1);

                            spanned_vertices.push_back(source(e,gm));
                            spanning_ids.push_back(spanning_ids.back());
                        }
                        else
                        {
                            spanning.push_back(e);
					        e_color[e] = GREEN;

                            spanned_vertices.push_back(target(e,gm));
                            spanning_ids.push_back(spanning_ids[pos1 - spanned_vertices.begin()]);
                        }
                    }
                    else if ( pos1 == spanned_vertices.end() )
                    {
                        spanning.push_back(e);
                        e_color[e] = GREEN;

                        spanned_vertices.push_back(source(e,gm));
                        spanning_ids.push_back(spanning_ids[pos2 - spanned_vertices.begin()]);
                    }
                    else
                    {
                        if ( spanning_ids[pos1 - spanned_vertices.begin()] < spanning_ids[pos2 - spanned_vertices.begin()] )
                        {
                            spanning.push_back(e);
                            e_color[e] = GREEN;

                            for (unsigned c=0, temp=spanning_ids[pos2 - spanned_vertices.begin()]; c < spanning_ids.size(); c++)
                                if (spanning_ids[c] == temp)
                                    spanning_ids[c] = spanning_ids[pos1 - spanned_vertices.begin()];
                        }
                        else if ( spanning_ids[pos1 - spanned_vertices.begin()] > spanning_ids[pos2 - spanned_vertices.begin()] )
                        {
                            spanning.push_back(e);
                            e_color[e] = GREEN;

                            for (unsigned c=0, temp=spanning_ids[pos1 - spanned_vertices.begin()]; c < spanning_ids.size(); c++)
                                if (spanning_ids[c] == temp)
                                    spanning_ids[c] = spanning_ids[pos2 - spanned_vertices.begin()];
                        }
                    }
                }
            }
        }
    }

    //fixer les pont
    {
        edge_iteratorM ei, ei_end;
        for(boost::tuples::tie(ei, ei_end) = edges(gm); ei != ei_end; ++ei)
        {
            edge_descriptorM e(*ei);

            if (e_color[e] == RED)
            {
                e_color[e] = GREEN;
                deque < vertex_descriptorM >::iterator pos1 = find(spanned_vertices.begin(), spanned_vertices.end(), source(e,gm));
                deque < vertex_descriptorM >::iterator pos2 = find(spanned_vertices.begin(), spanned_vertices.end(), target(e,gm));

                    if ( pos2 == spanned_vertices.end() )   //pos1 == spanned_vertices.end() ||
                    {
                        if ( pos1 == spanned_vertices.end() )
                        {
                            spanning.push_back(e);
					        e_color[e] = GREEN;

                            spanned_vertices.push_back(target(e,gm));
                            if(max_element(spanning_ids.begin(),spanning_ids.end()) == spanning_ids.end())
                                spanning_ids.push_back(0);
                            else
                                spanning_ids.push_back(*max_element(spanning_ids.begin(),spanning_ids.end())+1);

                            spanned_vertices.push_back(source(e,gm));
                            spanning_ids.push_back(spanning_ids.back());
                        }
                        else
                        {
                            spanning.push_back(e);
					        e_color[e] = GREEN;

                            spanned_vertices.push_back(target(e,gm));
                            spanning_ids.push_back(spanning_ids[pos1 - spanned_vertices.begin()]);
                        }
                    }
                    else if ( pos1 == spanned_vertices.end() )
                    {
                        spanning.push_back(e);
                        e_color[e] = GREEN;

                        spanned_vertices.push_back(source(e,gm));
                        spanning_ids.push_back(spanning_ids[pos2 - spanned_vertices.begin()]);
                    }
                    else
                    {
                        if ( spanning_ids[pos1 - spanned_vertices.begin()] < spanning_ids[pos2 - spanned_vertices.begin()] )
                        {
                            spanning.push_back(e);
                            e_color[e] = GREEN;

                            for (unsigned c=0, temp=spanning_ids[pos2 - spanned_vertices.begin()]; c < spanning_ids.size(); c++)
                                if (spanning_ids[c] == temp)
                                    spanning_ids[c] = spanning_ids[pos1 - spanned_vertices.begin()];
                        }
                        else if ( spanning_ids[pos1 - spanned_vertices.begin()] > spanning_ids[pos2 - spanned_vertices.begin()] )
                        {
                            spanning.push_back(e);
                            e_color[e] = GREEN;

                            for (unsigned c=0, temp=spanning_ids[pos1 - spanned_vertices.begin()]; c < spanning_ids.size(); c++)
                                if (spanning_ids[c] == temp)
                                    spanning_ids[c] = spanning_ids[pos2 - spanned_vertices.begin()];
                        }
                    }
            }
        }
    }

    //affichage des couleurs
    {
      /*  edge_iteratorM ei, ei_end;
        for(boost::tuples::tie(ei, ei_end) = edges(gm); ei != ei_end; ++ei)
        {
            cout << "arete d index : " << get(edge_index, gm, *ei)  << " entre les noeuds " << get(vertex_index, gm, source(*ei, gm)) << " et " << get(vertex_index, gm, target(*ei, gm));
            if(e_color[*ei] == BLACK)
                cout << " est de couleur BLACK" << endl;
            if(e_color[*ei] == GREY)
                cout << " est de couleur GREY" << endl;
            if(e_color[*ei] == GREEN)
                cout << " est de couleur GREEN" << endl;
            if(e_color[*ei] == RED)
                cout << " est de couleur RED" << endl;
        }
        cout << endl; */
    }

    //affichage des spanned vertices ainsi que leur arbre
    {
      /*  if(spanned_vertices.size()==spanning_ids.size())
            cout<<"taille identique  " << endl << endl;

        for(unsigned i(0); i < spanned_vertices.size();i++)
        {
            cout << "spanned vertex " << i << " d'index : " << get(vertex_index, gm, spanned_vertices[i])  << " appartient au bou de l'arbre " << spanning_ids[i] << endl;
        }   */
    }

    map < vertex_descriptorM, int > deg_r;

    //adjusting non corded cycles's vertices order
    {
        //map < vertex_descriptorM, int > deg_r;

        graph_traits<graphM>::out_edge_iterator out_i, out_end;

        for (unsigned i=0; i < g_vertices.size(); i++)
        {
            int deg(0);
            for (boost::tuples::tie(out_i, out_end) = out_edges(g_vertices[i], gm); out_i != out_end; ++out_i)
                if (e_color[*out_i] == BLACK)
                    deg++;
            deg_r[g_vertices[i]] = deg;
        }

        triFusionVerticesR(g_vertices_temp, g_vertices_temp.size(), gm, deg_r);
    }

    //affichage des vertices a traiter ainsi que leurs degrees relatifs
    {
       /* for (unsigned s(0) ; s < g_vertices_temp.size() ; s++)
            cout << "sommet " << s << " est d index : " << get(vertex_index, gm, g_vertices_temp[s])  << " est de degre relatif : " << deg_r[g_vertices_temp[s]]  << endl;
        cout << endl;  */
    }

    //calculating spanning for non corded cycles
    {
        graph_traits<graphM>::out_edge_iterator out_i, out_end;

        while ( g_vertices_temp.size() > 0 )
        {
            for (boost::tuples::tie(out_i, out_end) = out_edges(g_vertices_temp.front(), gm); out_i != out_end; ++out_i)
            {
                edge_descriptorM e(*out_i);
                e = edge(source(e,gm),target(e,gm),gm).first;

                if (e_color[e] == BLACK)
                {
                    deque < vertex_descriptorM >::iterator pos1 = find(spanned_vertices.begin(), spanned_vertices.end(), source(e,gm));    // or  g_vertices_temp.front());
                    deque < vertex_descriptorM >::iterator pos2 = find(spanned_vertices.begin(), spanned_vertices.end(), target(e,gm));

                    if ( pos2 == spanned_vertices.end() )   //pos1 == spanned_vertices.end() ||
                    {
                        if ( pos1 == spanned_vertices.end() )
				        {
					        spanning.push_back(e);
					        e_color[e] = GREEN;

                            spanned_vertices.push_back(target(e,gm));
					        spanning_ids.push_back(*max_element(spanning_ids.begin(),spanning_ids.end())+1);

                            spanned_vertices.push_back(source(e,gm));
					        spanning_ids.push_back(spanning_ids.back());
				        }
				        else
				        {
					        spanning.push_back(e);
					        e_color[e] = GREEN;

                            spanned_vertices.push_back(target(e,gm));
					        spanning_ids.push_back(spanning_ids[pos1 - spanned_vertices.begin()]);
				        }
                    }
                    else if ( pos1 == spanned_vertices.end() )
                    {
                        spanning.push_back(e);
                        e_color[e] = GREEN;

                        spanned_vertices.push_back(source(e,gm));
				        spanning_ids.push_back(spanning_ids[pos2 - spanned_vertices.begin()]);
                    }
			        else
			        {
				        if ( spanning_ids[pos1 - spanned_vertices.begin()] < spanning_ids[pos2 - spanned_vertices.begin()] )
				        {
					        spanning.push_back(e);
					        e_color[e] = GREEN;
					        for (unsigned c=0, temp=spanning_ids[pos2 - spanned_vertices.begin()]; c < spanning_ids.size(); c++)
						        if (spanning_ids[c] == temp)
							        spanning_ids[c] = spanning_ids[pos1 - spanned_vertices.begin()];
				        }
				        else if ( spanning_ids[pos1 - spanned_vertices.begin()] > spanning_ids[pos2 - spanned_vertices.begin()] )
				        {
					        spanning.push_back(e);
					        e_color[e] = GREEN;
					        for (unsigned c=0, temp=spanning_ids[pos1 - spanned_vertices.begin()]; c < spanning_ids.size(); c++)
						        if (spanning_ids[c] == temp)
							        spanning_ids[c] = spanning_ids[pos2 - spanned_vertices.begin()];
				        }
				        else
				        {
				            chords.push_back(e);
				            e_color[e] = GREY;
				        }
			        }
                }
            }
            g_vertices_temp.pop_front();
        }
    }

    //affichage des spanned vertices ainsi que leur arbre
    {
     /*   if(spanned_vertices.size()==spanning_ids.size())
            cout<<"taille identique  " << endl << endl;

        for(unsigned i(0); i < spanned_vertices.size();i++)
        {
            cout << "spanned vertex " << i << " d'index : " << get(vertex_index, gm, spanned_vertices[i])  << " appartient au bou de l'arbre " << spanning_ids[i] << endl;
        }  */
    }

    //affichage des aretes formant l'arbre de recouvrement ainsi que les cordes
    {
      /*  for (unsigned i =0 ; i < spanning.size() ; i++)
            cout << "Branche " << i << " est d index : " << get(edge_index, gm, spanning[i])  << " entre les noeuds " << get(vertex_index, gm, source(spanning[i], gm)) << " et " << get(vertex_index, gm, target(spanning[i], gm)) <<endl;
        cout << endl;

        for (unsigned i =0 ; i < chords.size() ; i++)
            cout << "Corde " << i << " est d index : " << get(edge_index, gm, chords[i])  << " entre les noeuds " << get(vertex_index, gm, source(chords[i], gm)) << " et " << get(vertex_index, gm, target(chords[i], gm)) <<endl;
        cout << endl; */
    }

}


//traitement d'initialisation des structures cycles et e_cycles
{
    for(unsigned i=0; i<cycles.size(); ++i)
        cycles[i].clear();

    for(unsigned i=0; i<e_cycles.size(); ++i)
        e_cycles[i].clear();

    cycles.clear();
    e_cycles.clear();

    for(unsigned i=0; i<chords.size(); ++i)
        cycles.push_back( deque < edge_descriptorM > () );

    for(unsigned i=0; i<num_edges(gm); ++i)
        e_cycles.push_back( deque <int> () );
}


//calcul des cycles
{

    tailleBCF=0;

    #pragma omp parallel for
    for (unsigned c=0; c < chords.size(); c++)
    {
        cycles[c].push_back(chords[c]);
        e_cycles[get(edge_index, gm, chords[c])].push_back(c);
        calculate_cycle(source(chords[c], gm), target(chords[c], gm), e_color, gm, cycles, e_cycles, c, e_index);

        tailleBCF += cycles[c].size();
    }

    cout << "La taille de la BCF est de " << tailleBCF << endl << endl;

    fichier_res << tailleBCF << endl;

}


    // get the end time
    gettimeofday(&endTime, NULL);
    // calculate time in microseconds
    tS = startTime.tv_sec*1000000 + (startTime.tv_usec);
    tE = endTime.tv_sec*1000000  + (endTime.tv_usec);
    std::cout << "Rendre la base de cycle fondamentale a dure : " << (tE - tS) / 1000000.0L << endl << endl;

    fichier_res << (tE - tS) / 1000000.0L << endl;


//fixing the gap between cycles and indexes
{
    for (unsigned i=0; i < cycles.size(); i++)
        for (unsigned j=0; j < cycles[i].size(); j++)
            cycles[i][j] = edge(source(cycles[i][j], gm), target(cycles[i][j], gm), gm).first;
}


//recalculer e_cycles en fonction de la nouvelle base de cycles
{
    for (unsigned i=0; i < e_cycles.size(); i++)
        e_cycles[i].clear();

    for (unsigned i=0; i < cycles.size(); i++)
        for (unsigned j=0; j < cycles[i].size(); j++)
            e_cycles[get(edge_index, gm, cycles[i][j])].push_back(i);
}


//affichage des cycles et en fonction de ces derniers, et en fonction des aretes
{
/*   for (unsigned i=0; i < cycles.size(); i++)
    {
        cout << "le cycle " << i << " est de taille " << cycles[i].size() << " et contient les aretes :" << endl;
        for (unsigned j=0; j < cycles[i].size(); j++)
        {
            cout << "    l'arete " << get(edge_index, gm, cycles[i][j]) << " reliant les sommets " << get(vertex_index, gm, source(cycles[i][j], gm)) << " et " << get(vertex_index, gm, target(cycles[i][j], gm)) << endl;
        }
    }
    cout << endl;


    for (unsigned i=0; i < e_cycles.size(); i++)
    {
        cout << "l'arete " << i << " appartient aux cycles :" << endl;
        for (unsigned j=0; j < e_cycles[i].size(); j++)
        {
            cout << "    le cycle " << e_cycles[i][j] << endl;
        }
    }
    cout << endl;   */

}


//verifier si la BC est fondamentale
{

    //edge_descriptorM e(get_edge(index, gm).first);

    for (unsigned i=0; i < e_cordes.size(); i++)
        e_cordes[i].clear();
    e_cordes.clear();


    //extraction des cordes probables
    for (unsigned i=0; i < e_cycles.size(); i++)
        if (e_cycles[i].size() == 1)
        {
            e_cordes.push_back(deque < int > ());
            e_cordes.back().push_back(i);
            e_cordes.back().push_back(e_cycles[i].front());
        }


    for (unsigned i=0; i < e_cordes.size(); i++)
    {
        unsigned j(i+1);
        while (j < e_cordes.size())
        {
            if (e_cordes[i][1] == e_cordes[j][1])
                e_cordes.erase(e_cordes.begin() + j);
            else
                j++;
        }
    }


    if (e_cordes.size() < cycles.size())
    {
        cout << " La BC n'est pas fondamentale !!! " << endl;
        cout << "Car " << cycles.size() - e_cordes.size() << " cycles ne possedent pas de cordes !" << endl;
    }
    else
        cout << "Le BC est FONDAMENTALE !!! " << endl;

}


//Calcul de la matrice d'adjacence des cycles
{

    // Create a 2D array
    typedef boost::multi_array<int, 2> matrix;
    //typedef matrix::index index;   //pour le parcours

    matrix mac(boost::extents[cycles.size()][cycles.size()]);

    for (unsigned i=0; i < e_cycles.size(); i++)
        for (unsigned j=0; j < e_cycles[i].size(); j++)
            for (unsigned k=j+1; k < e_cycles[i].size(); k++)
                mac[e_cycles[i][j]][e_cycles[i][k]] = 1;

    //Affichage de la MAC
    {
     /*   cout << endl;
        for (unsigned i=0; i < mac.size(); i++)
        {
            for (unsigned j=0; j < mac[i].size(); j++)
                cout << mac[i][j] << " ";
            cout << endl;
        } */
    }

    int numberOf1(0);

    //Calcul du nombre d'elements non nuls
    {
        for (unsigned i=0; i < mac.size(); i++)
            for (unsigned j=0; j < mac[i].size(); j++)
                if (mac[i][j]==1)
                    numberOf1++;
    }

    cout << endl << "Le nombre d'elements non nul est egal a " << numberOf1 << endl << endl;

    fichier_res << numberOf1 << endl;

    //verifier si triangulaire superieure
    {
        bool flag(true);
        for (unsigned i=0; i < mac.size(); i++)
            for (unsigned j=0; j < i; j++)
                if (mac[i][j]==1)
                    flag=false;
        if(flag)
            cout << "la matrice est triangulaire superieure." << endl << endl;
        else
            cout << "la matrice n'est pas triangulaire superieure." << endl << endl;
    }

}

    fichier_res.close();

    return 0;
}
