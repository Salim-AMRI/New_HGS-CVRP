#include "Genetic.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <fstream>
#include <unordered_map>

void Genetic::run()
{
        /* INITIAL POPULATION */
        population.generatePopulation();

        int cpt = 1;

        std::ofstream Evolutionfile("evolution/Instance_" + std::to_string(params.nbClients + 1) + "_" + "crossover_type_" + std::to_string(params.ap.crossoverType) + "_seed_" +  std::to_string(params.ap.seed) + ".txt");

        const double interval = 0.2; // Interval de temps en secondes
        double lastOutputTime = (double)clock() / (double)CLOCKS_PER_SEC; // Initialisation du temps de la dernière sortie

        int nbIter;
        int nbIterNonProd = 1;
        if (params.verbose) std::cout << "----- STARTING GENETIC ALGORITHM" << std::endl;

        std::cout << " params.ap.nbIter " <<   params.ap.nbIter << std::endl;
        std::cout << " params.ap.crossoverType " <<   params.ap.crossoverType << std::endl;

        for (nbIter = 0 ;  (double)(clock()-params.startTime)/(double)CLOCKS_PER_SEC < params.ap.timeLimit ; nbIter++)
        {
                /* SELECTION AND CROSSOVER */
                if (params.ap.crossoverType == 0) {
                    crossoverOX(offspring, population.getBinaryTournament(), population.getBinaryTournament());
                } else if (params.ap.crossoverType == 1) {
                    crossoverLOX(offspring, population.getBinaryTournament(), population.getBinaryTournament());
                } else if (params.ap.crossoverType == 2) {
                    crossoverAOX(offspring, population.getBinaryTournament(), population.getBinaryTournament());
                } else if (params.ap.crossoverType == 3) {
                    crossoverUOX(offspring, population.getBinaryTournament(), population.getBinaryTournament());
                } else if (params.ap.crossoverType == 4) {
                    crossoverPMX(offspring, population.getBinaryTournament(), population.getBinaryTournament());
                } else if (params.ap.crossoverType == 5) {
                    crossoverEAX(offspring, population.getBinaryTournament(), population.getBinaryTournament());
                } else if (params.ap.crossoverType == 6) {
                    crossoverGPX(offspring, population.getBinaryTournament(), population.getBinaryTournament());
                } else if (params.ap.crossoverType == 7) {
                    crossoverGrPX(offspring, population.getBinaryTournament(), population.getBinaryTournament());
                } else if (params.ap.crossoverType == 8) {
                    crossoverPathRelinking(offspring, population.getBinaryTournament(), population.getBinaryTournament());
                } else if (params.ap.crossoverType == 9) {
                    crossoverGrPX2(offspring, population.getBinaryTournament(), population.getBinaryTournament());
                } else if (params.ap.crossoverType == 10) {
                    crossoverASSC(offspring, population.getBinaryTournament(), population.getBinaryTournament());
                } else if (params.ap.crossoverType == 11) {
                    crossoverASSC2(offspring, population.getBinaryTournament(), population.getBinaryTournament());
                } else if (params.ap.crossoverType == 12) {
                    crossoverGOX(offspring, population.getBinaryTournament(), population.getBinaryTournament());
                }

        
        
        double currentTime = (double)clock() / (double)CLOCKS_PER_SEC;
        if (currentTime - lastOutputTime >= interval)
        {
            Evolutionfile << population.getBestFeasible()->eval.penalizedCost << " - " << currentTime << std::endl;
            lastOutputTime = currentTime;
        }


        /* LOCAL SEARCH */
        localSearch.run(offspring, params.penaltyCapacity, params.penaltyDuration);
                bool isNewBest = population.addIndividual(offspring,true);
                if (!offspring.eval.isFeasible && params.ran()%2 == 0) // Repair half of the solutions in case of infeasibility
                {
                        localSearch.run(offspring, params.penaltyCapacity*10., params.penaltyDuration*10.);
                        if (offspring.eval.isFeasible) isNewBest = (population.addIndividual(offspring,false) || isNewBest);
                }

                /* TRACKING THE NUMBER OF ITERATIONS SINCE LAST SOLUTION IMPROVEMENT */
                if (isNewBest) nbIterNonProd = 1;
                else nbIterNonProd ++ ;

                /* DIVERSIFICATION, PENALTY MANAGEMENT AND TRACES */
                if (nbIter % params.ap.nbIterPenaltyManagement == 0) population.managePenalties();
                if (nbIter % params.ap.nbIterTraces == 0) population.printState(nbIter, nbIterNonProd);

                /* FOR TESTS INVOLVING SUCCESSIVE RUNS UNTIL A TIME LIMIT: WE RESET THE ALGORITHM/POPULATION EACH TIME maxIterNonProd IS ATTAINED*/

// 		if (params.ap.timeLimit != 0 && nbIterNonProd == params.ap.nbIter)
// 		{
// 			population.restart();
// 			nbIterNonProd = 1;
//             MyFile << "Restart" << std::endl;
//
// 		}

        }


        Evolutionfile << population.getBestFeasible()->eval.penalizedCost << " - " << (double)(clock() - params.startTime) / (double)CLOCKS_PER_SEC << std::endl;

        if (params.verbose) std::cout << "----- GENETIC ALGORITHM FINISHED AFTER " << nbIter << " ITERATIONS. TIME SPENT: " << (double)(clock() - params.startTime) / (double)CLOCKS_PER_SEC << std::endl;

        Evolutionfile.close();
}

void Genetic::crossoverOX(Individual & result, const Individual & parent1, const Individual & parent2)
{
        // Frequency table to track the customers which have been already inserted
        std::vector <bool> freqClient = std::vector <bool> (params.nbClients + 1, false);

        // Picking the beginning and end of the crossover zone
        std::uniform_int_distribution<> distr(0, params.nbClients-1);
        int start = distr(params.ran);
        int end = distr(params.ran);

        // Avoid that start and end coincide by accident
        while (end == start) end = distr(params.ran);

        // Copy from start to end
        int j = start;
        while (j % params.nbClients != (end + 1) % params.nbClients)
        {
                result.chromT[j % params.nbClients] = parent1.chromT[j % params.nbClients];
                freqClient[result.chromT[j % params.nbClients]] = true;
                j++;
        }

        // Fill the remaining elements in the order given by the second parent
        for (int i = 1; i <= params.nbClients; i++)
        {
                int temp = parent2.chromT[(end + i) % params.nbClients];
                if (freqClient[temp] == false)
                {
                        result.chromT[j % params.nbClients] = temp;
                        j++;
                }
        }

        // Complete the individual with the Split algorithm
        split.generalSplit(result, parent1.eval.nbRoutes);
}

void Genetic::crossoverLOX(Individual & result, const Individual & parent1, const Individual & parent2)
{
    // Frequency table to track the customers which have been already inserted
    std::vector<bool> freqClient = std::vector<bool>(params.nbClients + 1, false);

    // Picking the beginning and end of the crossover zone
    std::uniform_int_distribution<> distr(0, params.nbClients - 1);
    int start = distr(params.ran);
    int end = distr(params.ran);

    // Ensure start is smaller than end
    if (start > end)
        std::swap(start, end);

    // Copy the crossover segment from parent1 to the offspring
    for (int i = start; i <= end; i++) {
        result.chromT[i % params.nbClients] = parent1.chromT[i % params.nbClients];
        freqClient[result.chromT[i % params.nbClients]] = true;
    }

    // Fill the remaining elements using parent2
    int j = (end + 1) % params.nbClients; // Start from the next index after end
    for (int i = 0; i < params.nbClients; i++) {
        if (!freqClient[parent2.chromT[i]]) {
            result.chromT[j] = parent2.chromT[i];
            freqClient[parent2.chromT[i]] = true;
            j = (j + 1) % params.nbClients;
        }
    }

    // Complete the individual with the Split algorithm
    split.generalSplit(result, parent1.eval.nbRoutes);
}

void Genetic::crossoverAOX(Individual & result, const Individual & p1, const Individual & p2)
{
    // Frequency table to track customers already inserted
    std::vector<bool> freqClient = std::vector<bool>(params.nbClients + 1, false);

    // Picking the beginning and end of the crossover zone
    int start1 = std::rand() % params.nbClients;
    int end1 = std::rand() % params.nbClients;
    while (end1 == start1) end1 = std::rand() % params.nbClients;

    // Shift zone in p2 to match final customer of zone in p1
    int start2 = start1, end2 = end1;
    while (p2.chromT[end2 % params.nbClients] != p1.chromT[end1 % params.nbClients])
        start2++, end2++;

    // Test if zone in p1 is different to zone in p2
    bool same = true;
    int size = (start1 < end1 ? end1 - start1 : params.nbClients - start1 + end1);
    for (int j = 0; j < size && same; j++) {
        if (p1.chromT[(start1 + j) % params.nbClients] != p2.chromT[(start2 + j) % params.nbClients])
            same = false;
    }

    // If same, randomize point in p2
    if (same) end2 = end2 + rand() % (params.nbClients - size);

    // Copy in place the elements from start to end
    int j = start1;
    while (j % params.nbClients != (end1 + 1) % params.nbClients) {
        result.chromT[j % params.nbClients] = p1.chromT[j % params.nbClients];
        freqClient[result.chromT[j % params.nbClients]] = true;
        j++;
    }

    // Fill the remaining elements in the order given by p2
    for (int i = 1; i <= params.nbClients; i++) {
        int temp = p2.chromT[(end2 + i) % params.nbClients];
        if (freqClient[temp] == false) {
            result.chromT[j % params.nbClients] = temp;
            j++;
        }
    }

    // Completing the individual with the Split algorithm
    split.generalSplit(result, p1.eval.nbRoutes);
}

void Genetic::crossoverUOX(Individual & result, const Individual & parent1, const Individual & parent2)
{
    // Frequency table to track the customers which have been already inserted
    std::vector<bool> freqClient = std::vector<bool>(params.nbClients + 1, false);

    // Picking crossover points
    std::uniform_int_distribution<> distr(0, params.nbClients - 1);
    int point1 = distr(params.ran);
    int point2 = distr(params.ran);

    // Ensure point1 and point2 are different
    while (point2 == point1)
        point2 = distr(params.ran);

    // Ensure point2 is greater than point1
    if (point1 > point2)
        std::swap(point1, point2);

    // Copy elements between points from parent1 to result
    for (int i = point1; i <= point2; ++i) {
        result.chromT[i] = parent1.chromT[i];
        freqClient[result.chromT[i]] = true;
    }

    // Fill remaining elements from parent2
    int j = (point2 + 1) % params.nbClients;
    for (int i = 0; i < params.nbClients; ++i) {
        int index = (point2 + 1 + i) % params.nbClients;
        int temp = parent2.chromT[index];
        if (!freqClient[temp]) {
            result.chromT[j] = temp;
            freqClient[temp] = true;
            j = (j + 1) % params.nbClients;
        }
    }

    // Complete the individual with the Split algorithm
    split.generalSplit(result, parent1.eval.nbRoutes);
}

void Genetic::crossoverPMX(Individual & result, const Individual & parent1, const Individual & parent2)
{
    // Picking the beginning and end of the crossover zone
    std::uniform_int_distribution<> distr(0, params.nbClients - 1);
    int start = distr(params.ran);
    int end = distr(params.ran);

    // Avoid that start and end coincide by accident
    while (end == start)
        end = distr(params.ran);

    // Copy the segment from parent1 to the result and create mapping table
    for (int i = start; i <= end; ++i) {
        result.chromT[i] = parent1.chromT[i];
    }

    // Create mapping table for elements in the segment
    std::unordered_map<int, int> mapping;
    for (int i = start; i <= end; ++i) {
        mapping[parent1.chromT[i]] = parent2.chromT[i];
    }

    // Map elements from parent2 to result, avoiding elements already in the segment
    for (int i = 0; i < params.nbClients; ++i) {
        if (i < start || i > end) {
            int value = parent2.chromT[i];
            while (mapping.find(value) != mapping.end()) {
                value = mapping[value];
            }
            result.chromT[i] = value;
        }
    }

    // Complete the individual with the Split algorithm
    split.generalSplit(result, parent1.eval.nbRoutes);
}

void Genetic::crossoverEAX(Individual & result, const Individual & parent1, const Individual & parent2)
{
    // Step 1: Create a AB consisting of all edges of parent1 and parent2
    std::vector<std::vector<int>> AB(params.nbClients, std::vector<int>(params.nbClients, 0));
    for (int i = 0; i < params.nbClients; ++i) {
        for (int j = 0; j < params.nbClients; ++j) {
            if (parent1.chromT[i] == j || parent2.chromT[i] == j) {
                AB[i][j] = 1;
            }
        }
    }

    // Step 2: Find closed cycles by tracing the edges of parent1 and parent2 alternately from GAB
    std::vector<bool> visited(params.nbClients, false);
    std::vector<std::vector<int>> cycles;
    for (int i = 0; i < params.nbClients; ++i) {
        if (!visited[i]) {
            std::vector<int> cycle;
            int currentNode = i;
            while (!visited[currentNode]) {
                cycle.push_back(currentNode);
                visited[currentNode] = true;
                int nextNode = -1;
                for (int j = 0; j < params.nbClients; ++j) {
                    if (!visited[j] && AB[currentNode][j] == 1) {
                        nextNode = j;
                        break;
                    }
                }
                if (nextNode == -1 || nextNode == i) {
                    break; // No unvisited neighbor found or cycle is complete
                }
                currentNode = nextNode;
            }
            if (!cycle.empty()) {
                cycles.push_back(cycle);
            }
        }
    }

    // Randomly shuffle the cycles
    std::random_shuffle(cycles.begin(), cycles.end());

    // Step 3: Construct Eset consisting of several AB-cycles selected randomly from all cycles
    std::vector<std::vector<int>> Eset;
    std::random_shuffle(cycles.begin(), cycles.end());
    for (const auto& cycle : cycles) {
        if (cycle.size() > 2) { // Exclure les cycles avec seulement deux arêtes
            // Vérifier et supprimer les clients ayant une valeur égale à 0
            std::vector<int> filteredCycle;
            for (int node : cycle) {
                if (node != 0) {
                    filteredCycle.push_back(node);
                }
            }
            if (!filteredCycle.empty()) {
                Eset.push_back(filteredCycle);
            }
        }
    }

    // Step 4: Manipulate edges of parent1 based on Eset
    for (const auto& cycle : Eset) {
        for (int i = 0; i < cycle.size(); ++i) {
            int currentNode = cycle[i];
            int nextNode = cycle[(i + 1) % cycle.size()];
            result.chromT[currentNode] = nextNode;
        }
    }

    // Step 5: Assemble sub-tours into a complete cycle
    bool completeCycle = false;
    while (!completeCycle) {
        int startNode = -1;
        for (int i = 0; i < params.nbClients; ++i) {
            if (!visited[i]) {
                startNode = i;
                break;
            }
        }
        if (startNode == -1) {
            completeCycle = true;
            break;
        }
        int currentNode = startNode;
        while (true) {
            visited[currentNode] = true;
            int nextNode = -1;
            for (int i = 0; i < params.nbClients; ++i) {
                if (!visited[i] && AB[currentNode][i] == 1) {
                    nextNode = i;
                    break;
                }
            }
            if (nextNode == -1) {
                result.chromT[currentNode] = startNode;
                break;
            }
            result.chromT[currentNode] = nextNode;
            currentNode = nextNode;
        }
    }

    // Déclaration d'un vecteur pour suivre les clients déjà vus dans le cycle final
    std::vector<bool> clientVus(params.nbClients, false);

    // Recherche des clients en double et leur remplacement
    for (int i = 0; i < params.nbClients; ++i) {
        int clientActuel = result.chromT[i];
        if (!clientVus[clientActuel]) {
            // Marquer le client comme vu
            clientVus[clientActuel] = true;
        } else {
            // Remplacer le client répété par le client manquant du parent1.chromT
            for (int j = 0; j < params.nbClients; ++j) {
                if (parent1.chromT[j] != 0 && !clientVus[parent1.chromT[j]]) {
                    result.chromT[i] = parent1.chromT[j];
                    clientVus[parent1.chromT[j]] = true;
                    break;
                }
            }
        }
    }

    // Complete the individual with the Split algorithm
    split.generalSplit(result, parent1.eval.nbRoutes);
}

void Genetic::crossoverGPX(Individual & result, const Individual & parent1, const Individual & parent2) {
    // Taille des chromosomes des parents
    int length = parent1.chromT.size();

    // Initialisation du chromosome de l'enfant résultant
    result.chromT.clear();

    // Définition de la matrice AB
    std::unordered_map<int, std::unordered_map<int, int>> AB;

    // Initialisation de la matrice AB
    for (int i = 0; i < length; ++i) {
        AB[i] = std::unordered_map<int, int>();
    }

    // Construction de la matrice AB qui combine les deux parents
    for (int i = 0; i < length - 1; ++i) {
        AB[parent1.chromT[i]][parent1.chromT[i + 1]] = 1; // Les arcs du parent 1
        AB[parent1.chromT[i + 1]][parent1.chromT[i]] = 1;

        // Vérifier si l'arête est partagée entre les deux parents
        int vertex1 = parent1.chromT[i];
        int vertex2 = parent1.chromT[i + 1];
        if (AB[parent2.chromT[i]][parent2.chromT[i + 1]] == 1 &&
            (parent2.chromT[i] == vertex1 || parent2.chromT[i] == vertex2) &&
            (parent2.chromT[i + 1] == vertex1 || parent2.chromT[i + 1] == vertex2)) {

            AB[vertex1][vertex2] = 2; // Les arcs partagés entre les deux parents
            AB[vertex2][vertex1] = 2;
        }
    }

    // Traitement pour le dernier arc du chromosome
    AB[parent1.chromT[length - 1]][parent1.chromT[0]] = 1;
    AB[parent1.chromT[0]][parent1.chromT[length - 1]] = 1;
    if (std::find(parent2.chromT.begin(), parent2.chromT.end(), parent1.chromT[length - 1]) != parent2.chromT.end() &&
        std::find(parent2.chromT.begin(), parent2.chromT.end(), parent1.chromT[0]) != parent2.chromT.end()) {
        AB[parent1.chromT[length - 1]][parent1.chromT[0]] = 2;
        AB[parent1.chromT[0]][parent1.chromT[length - 1]] = 2;
    }

    // Suppression des arêtes partagées (marquées avec la valeur 2) dans la matrice AB
    for (int i = 0; i < length; ++i) {
        for (int j = 0; j < length; ++j) {
            if (AB[i][j] == 2) {
                AB[i][j] = -1; // Marquer temporairement pour la suppression
                AB[j][i] = -1; // Marquer temporairement pour la suppression
            }
        }
    }

    // Suppression des arêtes temporairement marquées
    for (int i = 0; i < length; ++i) {
        for (int j = 0; j < length; ++j) {
            if (AB[i][j] == -1) {
                AB[i][j] = 0; // Suppression de l'arête
            }
        }
    }

    // Transformation de la matrice AB pour que tous les sommets aient un degré de 3
    for (int i = 0; i < length; ++i) {
        int degree = 0;
        for (int j = 0; j < length; ++j) {
            degree += AB[i][j]; // Calcul du degré du sommet i
        }
        for (int j = 0; j < length && degree < 3; ++j) {
            if (i != j && AB[i][j] == 0) {
                AB[i][j] = 1; // Ajouter une arête supplémentaire
                AB[j][i] = 1; // Arête bidirectionnelle
                degree++;
            }
        }
    }

    // Transformation de la matrice AB pour que tous les sommets aient un degré de 3
    int countDegreeThree = 0; // Initialisation du compteur
    for (int i = 0; i < length; ++i) {
        int degree = 0;
        for (int j = 0; j < length; ++j) {
            degree += AB[i][j]; // Calcul du degré du sommet i
        }
        if (degree == 3) {
            countDegreeThree++; // Incrémentation du compteur si le degré est égal à 3
        }
    }

    // Construction du chromosome de l'enfant résultant du crossover GPX
    std::vector<bool> visited(length, false); // Marque les sommets visités
    int current = 0; // Sommet de départ (peut être choisi arbitrairement)

    // Boucle de construction du chemin
    do {
        result.chromT.push_back(current); // Ajouter le sommet actuel à la séquence du chromosome de l'enfant
        visited[current] = true; // Marquer le sommet comme visité

        // Trouver le prochain sommet connecté non visité
        int next = -1;
        for (int i = 0; i < length; ++i) {
            if (AB[current][i] != 0 && !visited[i]) {
                next = i;
                break;
            }
        }

        // Si aucun sommet valide non visité n'est trouvé, choisir aléatoirement un sommet non visité et continuer la recherche des successeurs à partir de ce sommet
        if (next == -1) {
            for (int i = 0; i < length; ++i) {
                if (!visited[i]) {
                    next = i;
                    break;
                }
            }
        }

        // Si aucun sommet non visité n'est trouvé, revenir au sommet de départ
        if (next == -1) {
            next = 0; // Retourner au sommet de départ
        }

        current = next; // Définir le prochain sommet comme le sommet actuel

    } while (current != 0); // Continuer jusqu'à ce que nous revenions au sommet de départ

    // Suppression des sommets de valeur 0 dans le chromosome
    result.chromT.erase(std::remove(result.chromT.begin(), result.chromT.end(), 0), result.chromT.end());

    // Comparaison avec le parent 1 et ajout des clients manquants
    for (int i = 0; i < parent1.chromT.size(); ++i) {
        if (std::find(result.chromT.begin(), result.chromT.end(), parent1.chromT[i]) == result.chromT.end()) {
            result.chromT.push_back(parent1.chromT[i]); // Ajouter le client manquant au chromosome de l'enfant
        }
    }

    // Compléter l'individu avec l'algorithme Split
    split.generalSplit(result, parent1.eval.nbRoutes);
}

void Genetic::crossoverGrPX(Individual &result, const Individual &parent1, const Individual &parent2)
{
    int nbParents = 2;
    double routeCost[2][params.nbVehicles];
    double routeBenefit[2][params.nbVehicles];
    int tabAffectedClients[params.nbClients];

    for (int i=0; i<params.nbClients; i++){
        tabAffectedClients[i] = -1;
    }

    std::vector <Individual> vParents (2);

    vParents[0] = parent1;
    vParents[1] = parent2;


    for (int j=0; j<params.nbVehicles; j++){

        result.chromR[j].clear();
    }

    for (int k=0; k<nbParents; k++) {

        for (int l=0; l<params.nbVehicles; l++){
            routeCost[k][l]=0;
            routeBenefit[k][l]=0;
        }


        for (int l=0; l<params.nbVehicles; l++){


            int currentClient = 0;


            for (int customer : vParents[k].chromR[l]) {


                routeCost[k][l] += params.timeCost[currentClient][customer];
                routeBenefit[k][l] += params.timeCost[0][customer] * params.cli[customer].demand;



                currentClient = customer;




            }

            routeCost[k][l] += params.timeCost[0][currentClient];


        }


    }



    double valMax;
    int vehMax;

    int indice=0;

    for (int i=0; i<params.nbVehicles; i++) {



        Individual& currentParent = vParents[indice];
                double* currentRouteCost = routeCost[indice];
        double* currentRouteBenefit = routeBenefit[indice];

                valMax=-1;
                vehMax=(rand()/(double)RAND_MAX) * params.nbVehicles ;


        for (int j=0; j<params.nbVehicles; j++) {

            double currentVal;
            if(currentRouteCost[j] != 0){
                currentVal = (double)currentRouteBenefit[j]/ (double)currentRouteCost[j];
            } else {
                currentVal = -10;
            }


            if (currentVal>valMax) {
                valMax=currentVal;
                vehMax=j;
            }

        }



        for (int customer : currentParent.chromR[vehMax]) {


            if(tabAffectedClients[customer-1] < 0){

                tabAffectedClients[customer-1] = i;

                result.chromR[i].push_back(customer);

            }





        }


        routeCost[indice][vehMax]=0;
        routeBenefit[indice][vehMax]=0;


        int indice2;

        if(indice == 0){
            indice2 = 1;

        } else {
            indice2 = 0;
        }



        for (int l=0; l<params.nbVehicles; l++){
            routeCost[indice2][l]=0;
            routeBenefit[indice2][l]=0;
        }


        for (int l=0; l<params.nbVehicles; l++){


            int currentClient = 0;

            double totalDemande = 0;


            for (int customer : vParents[indice2].chromR[l]) {

                if(tabAffectedClients[customer-1] < 0){

                    routeCost[indice2][l] += params.timeCost[currentClient][customer] ;
                    routeBenefit[indice2][l] += params.timeCost[0][customer]* params.cli[customer].demand;

                    //routeBenefit[indice2][l] += 1;
                    //routeBenefit[indice2][l] += params.cli[customer].demand;

                    currentClient = customer;

                }


            }

            routeCost[indice2][l] += params.timeCost[0][currentClient];

        }

        if(indice == 0){
            indice = 1;

        } else {
            indice = 0;
        }
    }
    for (int i=0; i<params.nbClients; i++){

                if (tabAffectedClients[i]<0) {
                        int randVeh =(rand()/(double)RAND_MAX) * params.nbVehicles ;
            result.chromR[randVeh].push_back(i+1);
                }
        }

    result.evaluateCompleteCost(params);

}

std::vector<std::vector<int>> evalDistRoute(int r1, int r2, std::vector<int>  v1, std::vector<int>  v2){

    int c;

    int n = v1.size();
    int m = v2.size();


    std::vector<std::vector<int> > D (n+1, std::vector<int>(m+1));


    for (int i = 0; i < n+1; i++){
        D[i][0] = i;
    }


    for (int j = 0; j < m+1; j++){
        D[0][j] = j;
    }


    for (int i = 1; i < n+1; i++){

        for (int j = 1; j < m+1; j++){

            if (v1[i-1]  == v2[j-1]){

                c = 0;
            } else {

                c = 9999;

            }

            D[i][j] = std::min(std::min(D[i - 1][j] + 1, D[i][j - 1] + 1), D[i - 1][j - 1] + c);

        }

    }


    return D;
}


class Dmatrix {
    public:
        std::vector<std::vector<int> > tab;


};


void Genetic::crossoverPathRelinking(Individual & result, const Individual & parent1, const Individual & parent2)
{


        std::vector<std::vector<int> > A (params.nbVehicles, std::vector<int>(params.nbVehicles));

        std::vector <int> affectation = std::vector <int> (params.nbVehicles);
        std::vector <int> affectationBis = std::vector <int> (params.nbVehicles);

        std::vector<std::vector<int> > ifRoadReversed (params.nbVehicles, std::vector<int>(params.nbVehicles));



        std::vector<std::vector<Dmatrix>> storeDistMatrix(params.nbVehicles, std::vector<Dmatrix>(params.nbVehicles));




        for (int r1 = 0; r1 < params.nbVehicles; r1++){

            for (int r2 = 0; r2 < params.nbVehicles; r2++){


               std::vector<int> v1 = parent1.chromR[r1];
               std::vector<int> v2 = parent2.chromR[r2];

               int n = v1.size();
               int m = v2.size();


               std::vector<std::vector<int> > D1 = evalDistRoute(r1, r2, v1, v2);

               int d1 =  D1[n][m];

               std::reverse(v2.begin(),v2.end());



               std::vector<std::vector<int> > D2 = evalDistRoute(r1, r2, v1, v2);

               int d2 =  D2[n][m];



                if(d1 < d2){
                    A[r1][r2] = d1;
                    Dmatrix d;
                    d.tab = D1;

                    storeDistMatrix[r1][r2] = d;


                } else{
                    A[r1][r2] = d2;
                    ifRoadReversed[r1][r2] = 1;
                    Dmatrix d;
                    d.tab = D2;

                    storeDistMatrix[r1][r2] = d;
                }


            }

        }




        for (int c = 0; c < params.nbVehicles; c++){

            int minVal = 999999;
            int minI = -1;
            int minJ = -1;

            for (int i = 0; i < params.nbVehicles; i++){

                for (int j = 0; j < params.nbVehicles; j++){

                    if (A[i][j] < minVal){

                        minVal = A[i][j];
                        minI = i;
                        minJ = j;

                    }

                }

            }

            for (int k = 0; k < params.nbVehicles; k++){

                A[minI][k] = 9999;
                A[k][minJ] = 9999;

            }

            affectation[minI] = minJ;

            affectationBis[minJ] = minI;


        }


        for (int k = 0; k < params.nbVehicles; k++){
            int r1 = k;
            int r2 = int(affectation[k]);


        }


        std::vector<std::set <int>> toRemove (params.nbVehicles);

        std::vector<std::set <int>> toAdd (params.nbVehicles);

        std::vector<std::set <int>> F (params.nbVehicles);

        std::set <int> M;


        for (int k = 0; k < params.nbVehicles; k++){



            int cpt_move = 0;
            int r1 = k;
            int r2 = int(affectation[k]);

            std::vector<int> v1 = parent1.chromR[r1];
            std::vector<int> v2 = parent2.chromR[r2];

            int n = v1.size();
            int m = v2.size();



            if(ifRoadReversed[r1][r2] == 1){

//                 std::cout <<  "v2 Reversed " << std::endl ;

                std::reverse(v2.begin(),v2.end());

            }
            
            std::vector<std::vector<int> > D = storeDistMatrix[r1][r2].tab;


            int pos1 = n;
            int pos2 = m;


            int cpt = 0;




            while (pos1 != 0 || pos2!= 0){


                if(pos1 == 0){


                    toAdd[r2].insert(v2[pos2-1]);

                    M.insert(v2[pos2-1]);


                    pos2 = pos2 - 1;


                } else if(pos2 == 0) {

                    //M.push_back(v1[pos1-1]);

                    toRemove[r1].insert(v1[pos1-1]);

                    M.insert(v1[pos1-1]);

                   // M.push_back(0);

                    pos1 = pos1 - 1;


                } else if(D[pos1-1][pos2 - 1] == D[pos1][pos2]  && D[pos1-1][pos2 - 1] < D[pos1][pos2 - 1] && D[pos1-1][pos2-1] < D[pos1-1][pos2]){


                    F[r2].insert(v2[pos2-1]);

                    //F.push_back(0);
                    pos1 = pos1 - 1;
                    pos2 = pos2 - 1;




                } else if( D[pos1-1][pos2] < D[pos1-1][pos2 - 1]  && D[pos1-1][pos2] < D[pos1][pos2 - 1]){

                    //M.push_back(v1[pos1-1]);
                   // M.push_back(0);
                    toRemove[r1].insert(v1[pos1-1]);

                    M.insert(v1[pos1-1]);

                    pos1 = pos1 - 1;

                } else {

                    //M.push_back(v2[pos2-1]);

                    toAdd[r2].insert(v2[pos2-1]);

                    M.insert(v2[pos2-1]);


                    //M.push_back(0);
                    pos2 = pos2 - 1;

                }



            }



        }


        std::set< std::tuple<int,int,int>> relocateMoves;


        for (auto m : M) {
//             std::cout << "client " << m;

            for (int i = 0; i < params.nbVehicles; i++){

                if (toRemove[i].count(m)) {

                    for (int j = 0; j < params.nbVehicles; j++){

                        if (toAdd[j].count(m)) {

                            std::tuple <int, int, int> move = std::make_tuple(m, i, j);

                            relocateMoves.insert(move);
                        }

                    }

                }

            }

        }


        int nbRelocateMovesTotal = relocateMoves.size();
//         std::cout << "nbRelocateMovesTotal : " << nbRelocateMovesTotal << std::endl;

        int nbRelocateMovesTODO = nbRelocateMovesTotal/2;
//         std::cout << "nbRelocateMovesTODO : " << nbRelocateMovesTODO << std::endl;


        std::vector< std::tuple<int,int,int>> relocateMovesToPerform;

        std::sample(relocateMoves.begin(), relocateMoves.end(), std::back_inserter(relocateMovesToPerform), nbRelocateMovesTODO, std::mt19937{std::random_device{}()});


        std::random_shuffle(relocateMovesToPerform.begin(), relocateMovesToPerform.end());

        for (int k = 0; k < params.nbVehicles; k++){

            result.chromR[k] = parent1.chromR[k];

        }

        for (auto move : relocateMovesToPerform){

            int c = std::get<0>(move);
            int r1 = std::get<1>(move);
            int r2 = std::get<2>(move);



            result.chromR[r1].erase(std::remove(result.chromR[r1].begin(), result.chromR[r1].end(), c), result.chromR[r1].end());;



            std::vector<int> v2 = parent2.chromR[r2];

            if(ifRoadReversed[affectationBis[r2]][r2] == 1){

//                 std::cout <<  "v2 Reversed " << std::endl ;
                std::reverse(v2.begin(),v2.end());

            }


            int pivot = -1;

            bool notFound = true;

            int idx = 0;

            while(notFound){

                if (F[r2].count(v2[idx])) {

                    pivot = v2[idx];
                }

                if(v2[idx] == c){

                    notFound = false;
                }

                idx = idx+1;

            }

            F[r2].insert(c);

//             std::cout <<  "pivot " << pivot << std::endl ;


            notFound = true;
            idx = 0;

            if(pivot != -1){
                while(notFound){


                    if(result.chromR[affectationBis[r2]][idx] == pivot){

                        notFound = false;
                    }

                    idx = idx+1;

                }
            }

            result.chromR[affectationBis[r2]].insert(result.chromR[affectationBis[r2]].begin() + idx, c);



        }

        // Build up the rest of the Individual structure
        result.evaluateCompleteCost(params);

}

void Genetic::crossoverGrPX2(Individual &result, const Individual &parent1, const Individual &parent2)
{
    int nbParents = 2;
    double routeCost[2][params.nbVehicles];
    double routeBenefit[2][params.nbVehicles];
    int tabAffectedClients[params.nbClients];

    for (int i=0; i<params.nbClients; i++){
        tabAffectedClients[i] = -1;
    }

    std::vector <Individual> vParents (2);

    vParents[0] = parent1;
    vParents[1] = parent2;


    for (int j=0; j<params.nbVehicles; j++){

        result.chromR[j].clear();
    }

    for (int k=0; k<nbParents; k++) {

        for (int l=0; l<params.nbVehicles; l++){
            routeCost[k][l]=0;
            routeBenefit[k][l]=0;
        }


        for (int l=0; l<params.nbVehicles; l++){


            int currentClient = 0;


            for (int customer : vParents[k].chromR[l]) {


                routeCost[k][l] += params.timeCost[currentClient][customer];
                routeBenefit[k][l] += params.timeCost[0][customer] * params.cli[customer].demand;
                //routeBenefit[k][l] += params.timeCost[0][customer];

                currentClient = customer;




            }

            routeCost[k][l] += params.timeCost[0][currentClient];

        }


    }



    double valMax;
    int vehMax;

    int indice=0;

    int cpt = 0;

    for (int i=0; i<params.nbVehicles; i++) {



        Individual& currentParent = vParents[indice];
                double* currentRouteCost = routeCost[indice];
        double* currentRouteBenefit = routeBenefit[indice];

                valMax=-1;
                vehMax=(rand()/(double)RAND_MAX) * params.nbVehicles ;


        for (int j=0; j<params.nbVehicles; j++) {

            double currentVal;
            if(currentRouteCost[j] != 0){
                currentVal = (double)currentRouteBenefit[j]/ (double)currentRouteCost[j];
            } else {
                currentVal = -10;
            }


            if (currentVal>valMax) {
                valMax=currentVal;
                vehMax=j;
            }

        }



        for (int customer : currentParent.chromR[vehMax]) {


            if(tabAffectedClients[customer-1] < 0){

                tabAffectedClients[customer-1] = i;

                result.chromR[i].push_back(customer);

            }





        }


        routeCost[indice][vehMax]=0;
        routeBenefit[indice][vehMax]=0;


        int indice2;

        if(indice == 0){
            indice2 = 1;

        } else {
            indice2 = 0;
        }



        for (int l=0; l<params.nbVehicles; l++){
            routeCost[indice2][l]=0;
            routeBenefit[indice2][l]=0;
        }


        for (int l=0; l<params.nbVehicles; l++){


            int currentClient = 0;


            for (int customer : vParents[indice2].chromR[l]) {

                if(tabAffectedClients[customer-1] < 0){

                    routeCost[indice2][l] += params.timeCost[currentClient][customer];
                    routeBenefit[indice2][l] += params.timeCost[0][customer]* params.cli[customer].demand;
                    //routeBenefit[indice2][l] += params.timeCost[0][customer];
                    currentClient = customer;

                }


            }

            routeCost[indice2][l] += params.timeCost[0][currentClient];

        }

        if(indice == 0 and cpt == 1){
            indice = 1;

        } else if(indice == 0 and cpt == 0){
            cpt = 1;

        } else if(indice == 1){
            indice = 0;
            cpt = 0;
        }
    }


    for (int i=0; i<params.nbClients; i++){

                if (tabAffectedClients[i]<0) {
                        int randVeh =(rand()/(double)RAND_MAX) * params.nbVehicles ;
            result.chromR[randVeh].push_back(i+1);
                }
        }

    result.evaluateCompleteCost(params);

}

void Genetic::crossoverASSC(Individual &result, const Individual &parent1, const Individual &parent2) {
    // Vecteur pour suivre la fréquence des clients déjà sélectionnés
    std::vector<bool> freqClient(params.nbClients + 1, false);

    // Point de coupe
    int pointCoupe = params.ap.nbCut;
    

    
    // Variables booléennes pour les options de crossover
    bool segmentsEqual = false;
    if(params.ap.eqSeg == 1){
        segmentsEqual = true;
    }
    
    bool useCostBenefit = false;
    if(params.ap.useCostBenefit == 1){
        useCostBenefit = true;
    }

    bool randomSegmentSelection = false;
    if(params.ap.randSelect == 1){
        randomSegmentSelection = true;
    }
    
    bool insertionSegment = false;
    if(params.ap.insertSeg == 1){
        insertionSegment = true;
    }
    
    

    // Génération aléatoire des points de coupe
    std::vector<int> pointsDeCoupeAleatoires;
    std::uniform_int_distribution<int> distr(1, params.nbClients - 2); // Exclut le client 0 et le client nbClients - 1
    std::mt19937 gen(std::random_device{}());

    if (segmentsEqual) {
        // Segments de taille uniforme
        while (pointsDeCoupeAleatoires.size() < pointCoupe) {
            int point = distr(gen);
            if (std::find(pointsDeCoupeAleatoires.begin(), pointsDeCoupeAleatoires.end(), point) == pointsDeCoupeAleatoires.end()) {
                pointsDeCoupeAleatoires.push_back(point);
            }
        }
    } else {  
        // Génération des points de coupe non uniforme
        int firstPoint = distr(gen); // Premier point de coupe aléatoire
        pointsDeCoupeAleatoires.push_back(firstPoint);

        // Calcul de la taille des segments
        int segmentSize = (params.nbClients - 1) / pointCoupe;

        // Génération des points de coupe suivants
        for (int i = 1; i < pointCoupe; ++i) {
            int nextPoint = pointsDeCoupeAleatoires.back() + segmentSize + 1;

            // Réajustement si le point dépasse les limites
            while (nextPoint >= params.nbClients) {
                nextPoint -= params.nbClients - 1;
            }

            pointsDeCoupeAleatoires.push_back(nextPoint);
        }
    }

    // Tri des points de coupe
    std::sort(pointsDeCoupeAleatoires.begin(), pointsDeCoupeAleatoires.end());

    // Vecteur de vecteurs pour stocker les segments
    std::vector<std::vector<int>> segments(pointsDeCoupeAleatoires.size());

    // Initialisation des compteurs
    size_t cpt1 = 0; // Indice du point de coupe
    size_t cpt2 = pointsDeCoupeAleatoires[0]; // Indice du client dans chromT
    size_t cpt3 = 0; // Indice du segment

    // Vecteur pour garder une trace des clients déjà inclus dans les segments précédents
    std::vector<bool> clientDejaAjoute(params.nbClients + 1, false);

    // Boucle de remplissage des segments
    while (cpt3 < segments.size()) {
        int currentClient = parent1.chromT[cpt2];

        // Ajouter le client au segment courant s'il n'est pas déjà ajouté
        if (!clientDejaAjoute[currentClient]) {
            segments[cpt3].push_back(currentClient);
            clientDejaAjoute[currentClient] = true;
        }

        cpt2 = (cpt2 + 1) % parent1.chromT.size(); // Passage au client suivant ou retour au début

        // Vérifier si on a atteint le prochain point de coupe
        if (cpt1 + 1 < pointsDeCoupeAleatoires.size() && cpt2 == pointsDeCoupeAleatoires[cpt1 + 1]) {
            // Passage au prochain segment si un point de coupe est atteint
            cpt1++;
            cpt3++;
        }

        // Si on revient au premier point de coupe après avoir atteint la fin, arrêter le remplissage
        if (cpt3 == segments.size() - 1 && cpt2 == pointsDeCoupeAleatoires[0]) {
            break;
        }
    }

    // Fonction locale pour évaluer un segment en fonction des critères choisis
    auto evaluateSegment = [&](const std::vector<int>& segment, const Individual& parent) -> double {
        double totalCost = 0.0;
        double totalBenefit = 0.0;

        for (const auto& parent : {parent1, parent2}) {
            for (int l = 0; l < params.nbVehicles; ++l) {
                double routeCost = 0.0;
                double routeBenefit = 0.0;
                int currentClient = 0;

                for (int customer : segment) {
                    if (!freqClient[customer]) {
                        routeCost += params.timeCost[currentClient][customer];
                        routeBenefit += params.timeCost[0][customer] * params.cli[customer].demand;
                        currentClient = customer;
                        freqClient[customer] = true;
                    }
                }

                totalCost += routeCost + params.timeCost[currentClient][0];
                totalBenefit += routeBenefit;
            }
        }

        return useCostBenefit ? totalBenefit : totalCost;
    };

    // Évaluation des segments générés
    std::vector<double> scores1(pointCoupe), scores2(pointCoupe);
    for (int i = 0; i < pointCoupe; ++i) {
        scores1[i] = evaluateSegment(segments[i], parent1);
        scores2[i] = evaluateSegment(segments[i], parent2);
    }

    // Liste pour stocker les segments sélectionnés
    std::vector<std::vector<int>> selectedSegmentsList;
    std::vector<bool> segmentSelected(pointCoupe, false); // Pour suivre les segments déjà sélectionnés
    std::unordered_set<int> newlySelectedClients; // Pour suivre les clients nouvellement sélectionnés

    // Sélection des segments en fonction des scores
    for (int i = 0; i < pointCoupe; ++i) {
        int parentIndex = (i % 2 == 0) ? 0 : 1; // Alternance entre les parents

        int minScoreSegmentIndex = -1;
        double minScore = std::numeric_limits<double>::max();

        if (randomSegmentSelection) {
            // Sélection aléatoire d'un segment non sélectionné
            do {
                minScoreSegmentIndex = std::rand() % pointCoupe;
            } while (segmentSelected[minScoreSegmentIndex]);
        } else {
            // Sélection déterministe du segment avec le score minimum
            for (int j = 0; j < pointCoupe; ++j) {
                if (!segmentSelected[j] && (parentIndex == 0 ? scores1[j] : scores2[j]) < minScore) {
                    minScoreSegmentIndex = j;
                    minScore = (parentIndex == 0 ? scores1[j] : scores2[j]);
                }
            }
        }

        // Si aucun segment n'est trouvé, sortir de la boucle
        if (minScoreSegmentIndex == -1) {
            break;
        }

        // Sélection du segment avec le score minimum
        selectedSegmentsList.push_back(segments[minScoreSegmentIndex]);

        // Mise à jour de freqClient pour marquer les clients du segment comme déjà sélectionnés
        for (int client : segments[minScoreSegmentIndex]) {
            freqClient[client] = true;
            newlySelectedClients.insert(client);
        }

        // Marquer le segment comme sélectionné
        segmentSelected[minScoreSegmentIndex] = true;

        // Mettre à jour les scores des segments restants affectés par la nouvelle sélection
        for (int j = 0; j < pointCoupe; ++j) {
            if (!segmentSelected[j]) {
                bool affected = false;
                for (int client : segments[j]) {
                    if (newlySelectedClients.find(client) != newlySelectedClients.end()) {
                        affected = true;
                        break;
                    }
                }

                if (affected) {
                    scores1[j] = evaluateSegment(segments[j], parent1);
                    scores2[j] = evaluateSegment(segments[j], parent2);
                }
            }
        }
    }
    
    // Mélanger les segments et les insérer de manière aléatoire
    if (!insertionSegment) {
        std::vector<int> segmentIndices(selectedSegmentsList.size());
        std::iota(segmentIndices.begin(), segmentIndices.end(), 0);
        std::shuffle(segmentIndices.begin(), segmentIndices.end(), std::mt19937{std::random_device{}()});

        for (size_t i = 0; i < segmentIndices.size(); ++i) {
            int segmentIndex = segmentIndices[i];
            // Insérer les clients du segment sélectionné dans l'offspring
            for (size_t j = 0; j < selectedSegmentsList[segmentIndex].size(); ++j) {
                int client = selectedSegmentsList[segmentIndex][j];
                if (!freqClient[client]) {
                    result.chromT.push_back(client);
                    freqClient[client] = true;
                }
            }
        }
    } else {
        // Si l'insertion de segment est activée, insérez les segments de manière déterministe
        std::vector<bool> segmentSelected(segments.size(), false);

        int lastInsertedClient = 0; // Initialiser le dernier client inséré à l'entrepôt
        for (size_t i = 0; i < selectedSegmentsList.size(); ++i) {
            // Trouver le segment avec le premier client le plus proche du dernier client inséré
            double minDistance = std::numeric_limits<double>::max();
            size_t nextSegmentIndex = 0;
            for (size_t j = 0; j < selectedSegmentsList.size(); ++j) {
                if (!segmentSelected[j]) { // Vérifier si le segment n'a pas déjà été sélectionné
                    double distance = params.timeCost[lastInsertedClient][parent1.chromT[selectedSegmentsList[j].front()]];
                    if (distance < minDistance) {
                        minDistance = distance;
                        nextSegmentIndex = j;
                    }
                }
            }

            // Insérer les clients du segment sélectionné dans l'offspring
            for (int client : segments[nextSegmentIndex]) {
                int selectedParentIndex = (selectedSegmentsList.size() % 2 == 0) ? 0 : 1;
                int clientIndex = (selectedParentIndex == 0) ? parent1.chromT[client] : parent2.chromT[client];
                if (!freqClient[clientIndex]) {
                    result.chromT.push_back(clientIndex);
                    freqClient[clientIndex] = true;
                }
            }

            // Marquer le segment comme sélectionné
            segmentSelected[nextSegmentIndex] = true;
        }
    }
   
    // Compléter l'individu avec l'algorithme de Split
    split.generalSplit(result, parent1.eval.nbRoutes);
}

void Genetic::crossoverASSC2(Individual &result, const Individual &parent1, const Individual &parent2) {
    // Vecteur pour suivre la fréquence des clients déjà sélectionnés
    std::vector<bool> freqClient(params.nbClients + 1, false);

    // Point de coupe initialisé à 2
    int pointCoupe = 3;
    // Nombre de segments à créer
    int numSegment = pointCoupe + 1;
    // Variable Booleen
    bool segmentsEqual = true;
    bool useCostBenefit = true;
    bool randomSegmentSelection = true;
    bool insertionSegment = true;

    // Vecteur pour stocker les segments
    std::vector<std::pair<int, int>> segments;

    // Générer les segments selon la valeur de segmentsEqual
    if (!segmentsEqual) {
        // Segments de taille égale
        int segmentLength = params.nbClients / numSegment;
        int remainder = params.nbClients % numSegment;
        int start = 0;

        for (int i = 0; i < numSegment; ++i) {
            int end = start + segmentLength - 1;
            if (i < remainder) { // Distribuer le reste parmi les premiers segments
                end++;
            }
            segments.push_back({start, end});
            start = end + 1;
        }
    } else {
        // Segments de taille variable
        std::vector<int> cutPoints(numSegment - 1);
        std::iota(cutPoints.begin(), cutPoints.end(), 1); // [1, 2, ..., numSegment-1]
        std::shuffle(cutPoints.begin(), cutPoints.end(), std::mt19937{std::random_device{}()});
        cutPoints.resize(numSegment - 1);

        // Ajouter les points de départ (0) et de fin (params.nbClients)
        cutPoints.insert(cutPoints.begin(), 0);
        cutPoints.push_back(params.nbClients);

        // Trier les points de coupe pour former les segments
        std::sort(cutPoints.begin(), cutPoints.end());

        for (int i = 0; i < numSegment; ++i) {
            segments.push_back({cutPoints[i], cutPoints[i + 1] - 1});
        }
    }
    
    int nbParents = 2;

    std::vector <Individual> vParents (2);

    vParents[0] = parent1;
    vParents[1] = parent2;

    // Définir la fonction pour évaluer les segments en tenant compte du bénéfice ou non
    auto evaluateSegment = [&](const std::pair<int, int>& segment, const std::vector<int>& chrom, const std::vector<bool>& freqClient) -> double {
        if (!useCostBenefit) {
            // Évaluation par coût uniquement
            double totalDistance = 0.0;
            for (int i = segment.first; i < segment.second; ++i) {
                int clientA = chrom[i];
                int clientB = chrom[i + 1];
                totalDistance += params.timeCost[clientA][clientB];
            }
            return totalDistance;
        } else {
            // Évaluation par coût bénéfice
            double totalCost = 0.0;
            double totalBenefit = 0.0;

            // Calculer le coût et le bénéfice pour chaque parent
            for (int k = 0; k < nbParents; ++k) {
                for (int l = 0; l < params.nbVehicles; ++l) {
                    double routeCost = 0.0;
                    double routeBenefit = 0.0;
                    int currentClient = 0;

                    for (int i = segment.first; i <= segment.second; ++i) {
                        int customer = chrom[i];
                        routeCost += params.timeCost[currentClient][customer];
                        routeBenefit += params.timeCost[0][customer] * params.cli[customer].demand;
                        currentClient = customer;
                    }

                    // Ajouter le coût et le bénéfice de la route à l'ensemble
                    totalCost += routeCost;
                    totalBenefit += routeBenefit;
                    
                    // Ajouter le coût de retour à l'entrepôt
                    totalCost += params.timeCost[currentClient][0];
                }
            }

            // Vous pouvez choisir comment combiner le coût et le bénéfice ici
            return totalCost;
            //return totalCost - totalBenefit;
        }
    };

    // Vecteurs pour stocker les scores de chaque segment pour les deux parents
    std::vector<double> scores1(numSegment), scores2(numSegment);


    for (int i = 0; i < numSegment; i++) {
        scores1[i] = evaluateSegment(segments[i], parent1.chromT, freqClient);
        scores2[i] = evaluateSegment(segments[i], parent2.chromT, freqClient);
    }

    // Vecteur pour stocker les segments sélectionnés et les parents utilisés
    std::vector<int> selectedSegments(numSegment, -1);
    std::vector<bool> parentUsed(2, false); // Pour suivre quels parents ont été utilisés
    
    // Sélectionner les segments alternativement des parents et les stocker dans une liste
    std::vector<std::pair<int, int>> selectedSegmentsList;
    std::mt19937 rng{std::random_device{}()}; // Générateur de nombres aléatoires
    for (int i = 0; i < numSegment; ++i) {
        int bestSegment = -1;
        double minDistance = std::numeric_limits<double>::max();
        int parentIndex = -1;
        
        if (!randomSegmentSelection) {
            // Sélection aléatoire d'un segment non sélectionné
            std::vector<int> availableSegments;
            for (int j = 0; j < numSegment; ++j) {
                if (selectedSegments[j] == -1) {
                    availableSegments.push_back(j);
                }
            }
            std::uniform_int_distribution<int> dist(0, availableSegments.size() - 1);
            bestSegment = availableSegments[dist(rng)];
            parentIndex = (i % 2 == 0) ? 0 : 1;
        } else {
            // Sélection déterministe du segment avec le score minimum
            for (int j = 0; j < numSegment; ++j) {
                if (selectedSegments[j] == -1) {
                    double score = (i % 2 == 0) ? scores1[j] : scores2[j];
                    if (bestSegment == -1 || score < minDistance) {
                        bestSegment = j;
                        minDistance = score;
                        parentIndex = (i % 2 == 0) ? 0 : 1;
                    }
                }
            }
        }
        
        selectedSegments[bestSegment] = parentIndex;
        parentUsed[parentIndex] = true;
        selectedSegmentsList.push_back(segments[bestSegment]);
    }
    
    if(!insertionSegment){
        // Si l'insertion de segment est activée, mélangez les segments sélectionnés
        // Mélanger les segments et les insérer de manière aléatoire
        std::vector<int> segmentIndices(segments.size());
        std::iota(segmentIndices.begin(), segmentIndices.end(), 0);
        std::shuffle(segmentIndices.begin(), segmentIndices.end(), std::mt19937{std::random_device{}()});

        for (size_t i = 0; i < segmentIndices.size(); ++i) {
            int segmentIndex = segmentIndices[i];
            // Insérer les clients du segment sélectionné dans l'offspring
            for (int j = segments[segmentIndex].first; j <= segments[segmentIndex].second; ++j) {
                int client = (selectedSegments[segmentIndex] == 0) ? parent1.chromT[j] : parent2.chromT[j];
                if (!freqClient[client]) {
                    result.chromT.push_back(client);
                    freqClient[client] = true;
                }
            }
        }
    } else {
    
        // Initialiser la liste des segments sélectionnés
        std::vector<bool> segmentSelected(segments.size(), false);

        // Insérer les segments dans l'offspring
        int lastInsertedClient = 0; // initialiser le dernier client inséré à l'entrepôt
        for (size_t i = 0; i < segments.size(); ++i) {
            // Trouver le segment avec le premier client le plus proche du dernier client inséré
            double minDistance = std::numeric_limits<double>::max();
            size_t nextSegmentIndex = 0;
            for (size_t j = 0; j < segments.size(); ++j) {
                if (!segmentSelected[j]) { // vérifier si le segment n'a pas été déjà sélectionné
                    double distance = params.timeCost[lastInsertedClient][parent1.chromT[segments[j].first]];
                    if (distance < minDistance) {
                        minDistance = distance;
                        nextSegmentIndex = j;
                    }
                }
            }

            // Insérer les clients du segment sélectionné dans l'offspring
            for (int j = segments[nextSegmentIndex].first; j <= segments[nextSegmentIndex].second; ++j) {
                int client = parent1.chromT[j];
                if (!freqClient[client]) {
                    result.chromT[j] = client;
                    freqClient[client] = true;
                    lastInsertedClient = client; // mettre à jour le dernier client inséré
                }
            }

            // Marquer le segment comme sélectionné
            segmentSelected[nextSegmentIndex] = true;
        }
    }
    
    // Compléter l'individu avec l'algorithme de Split
    split.generalSplit(result, parent1.eval.nbRoutes);
}

void Genetic::crossoverGOX(Individual &result, const Individual &parent1, const Individual &parent2)
{
    int nbClients = params.nbClients;
    
    // Point de coupe
    int pointCut = params.ap.nbCut;
    
    // Variables booléennes pour les options de crossover
    bool segmentsEqual = false;
    if(params.ap.eqSeg == 1){
        segmentsEqual = true;
    }
    
    bool useCostBenefit = false;
    if(params.ap.useCostBenefit == 1){
        useCostBenefit = true;
    }

    bool randomSegmentSelection = false;
    if(params.ap.randSelect == 1){
        randomSegmentSelection = true;
    }
    
    int insertionSegment = params.ap.insertSeg;

    // Vérification que pointCut est valide
    if (pointCut < 2 || pointCut > nbClients) {
        throw std::invalid_argument("pointCut doit être entre 2 et nbClients.");
    }

    std::vector<int> cutPoints(pointCut);

    if (segmentsEqual) {
        std::uniform_int_distribution<> distr(0, nbClients - 1);
        std::set<int> cuts;
        while (cuts.size() < pointCut) {
            cuts.insert(distr(params.ran));
        }
        cutPoints.assign(cuts.begin(), cuts.end());
    } else {
        int step = nbClients / pointCut;
        for (int i = 0; i < pointCut; ++i) {
            cutPoints[i] = i * step;
        }
    }

    std::sort(cutPoints.begin(), cutPoints.end());

    // Initialiser un tableau pour suivre les éléments déjà ajoutés
    std::vector<bool> freqClient(nbClients + 1, false);

    // Lambda pour évaluer un segment
    auto evaluateSegment = [&](const std::vector<int>& segment, const Individual& parent) -> double {
        double totalCost = 0.0;
        double totalBenefit = 0.0;
        std::vector<bool> tempFreqClient(freqClient);

        int currentClient = 0;
        for (int customer : segment) {
            if (!tempFreqClient[customer]) {
                double cost = params.timeCost[currentClient][customer];
                double benefit = params.timeCost[0][customer] * params.cli[customer].demand;
                totalCost += cost;
                totalBenefit += benefit;
                currentClient = customer;
                tempFreqClient[customer] = true;
            }
        }
        totalCost += params.timeCost[currentClient][0];
        return useCostBenefit ? totalBenefit : totalCost;
    };

    // Création et évaluation des segments pour les deux parents
    std::vector<std::vector<int>> segments(pointCut);
    for (int k = 0; k < pointCut; ++k) {
        int start = cutPoints[k];
        int end = (k + 1 < pointCut) ? cutPoints[k + 1] : nbClients;
        for (int i = start; i < end; ++i) {
            segments[k].push_back(parent1.chromT[i % nbClients]);
        }
    }

    std::vector<double> scores1(pointCut), scores2(pointCut);
    for (int i = 0; i < pointCut; ++i) {
        scores1[i] = evaluateSegment(segments[i], parent1);
        scores2[i] = evaluateSegment(segments[i], parent2);
    }

    // Sélection des segments
    std::vector<std::vector<int>> selectedSegmentsList;
    std::vector<bool> segmentSelected(pointCut, false);
    std::unordered_set<int> newlySelectedClients;

    for (int i = 0; i < pointCut; ++i) {
        int parentIndex = (i % 2 == 0) ? 0 : 1;
        int minScoreSegmentIndex = -1;
        double minScore = std::numeric_limits<double>::max();

        if (randomSegmentSelection) {
            do {
                minScoreSegmentIndex = std::rand() % pointCut;
            } while (segmentSelected[minScoreSegmentIndex]);
        } else {
            for (int j = 0; j < pointCut; ++j) {
                double score = (parentIndex == 0 ? scores1[j] : scores2[j]);
                if (!segmentSelected[j] && score < minScore) {
                    minScoreSegmentIndex = j;
                    minScore = score;
                }
            }
        }

        if (minScoreSegmentIndex == -1) break;

        selectedSegmentsList.push_back(segments[minScoreSegmentIndex]);
        for (int client : segments[minScoreSegmentIndex]) {
            freqClient[client] = true;
            newlySelectedClients.insert(client);
        }
        segmentSelected[minScoreSegmentIndex] = true;

        for (int j = 0; j < pointCut; ++j) {
            if (!segmentSelected[j]) {
                bool affected = false;
                for (int client : segments[j]) {
                    if (newlySelectedClients.find(client) != newlySelectedClients.end()) {
                        affected = true;
                        break;
                    }
                }
                if (affected) {
                    scores1[j] = evaluateSegment(segments[j], parent1);
                    scores2[j] = evaluateSegment(segments[j], parent2);
                }
            }
        }
    }

    // Gestion de l'insertion des segments dans le résultat
    if (insertionSegment == 0) {
        // Insertion simple des segments
        int j = 0;
        for (const auto& segment : selectedSegmentsList) {
            for (int client : segment) {
                result.chromT[j++] = client;
            }
        }
        // Compléter avec les clients restants du parent2
        for (int i = 0; i < nbClients; ++i) {
            int client = parent2.chromT[i];
            if (!freqClient[client]) {
                result.chromT[j++] = client;
                freqClient[client] = true;
            }
        }
    } else if (insertionSegment == 1) {
        // Insertion aléatoire des segments
        std::shuffle(selectedSegmentsList.begin(), selectedSegmentsList.end(), std::mt19937{std::random_device{}()});
        int j = 0;
        for (const auto& segment : selectedSegmentsList) {
            for (int client : segment) {
                result.chromT[j++] = client;
            }
        }
        // Compléter avec les clients restants du parent2
        for (int i = 0; i < nbClients; ++i) {
            int client = parent2.chromT[i];
            if (!freqClient[client]) {
                result.chromT[j++] = client;
            }
        }
    } else if (insertionSegment == 2) {
        
        // Méthode 2: Insertion avancée optimisée
        // Pré-calculer les distances entre les clients
        std::vector<std::vector<double>> distances(nbClients + 1, std::vector<double>(nbClients + 1, 0.0));
        for (int i = 0; i <= nbClients; ++i) {
            for (int j = 0; j <= nbClients; ++j) {
                distances[i][j] = params.timeCost[i][j];
            }
        }

        // Créer un tableau pour stocker les clients déjà ajoutés
        std::vector<bool> added(nbClients + 1, false);

        // Initialiser le dernier client inséré
        int lastInsertedClient = 0;
        result.chromT.clear();  // Réinitialiser le chromosome résultant

        // Priorité des segments à insérer en fonction des distances minimales
        std::vector<size_t> segmentOrder(selectedSegmentsList.size());
        std::iota(segmentOrder.begin(), segmentOrder.end(), 0); // Remplir avec 0, 1, 2, ..., numCuts - 1
        std::sort(segmentOrder.begin(), segmentOrder.end(), [&](size_t a, size_t b) {
            double distA = distances[lastInsertedClient][selectedSegmentsList[a].front()];
            double distB = distances[lastInsertedClient][selectedSegmentsList[b].front()];
            return distA < distB;
        });

        // Insérer les segments sélectionnés
        for (size_t i : segmentOrder) {
            const auto& segment = selectedSegmentsList[i];
            for (int client : segment) {
                if (!added[client]) {
                    result.chromT.push_back(client);
                    added[client] = true;
                    lastInsertedClient = client;
                }
            }
        }

        // Compléter avec les clients restants du parent2
        for (int i = 0; i < nbClients; ++i) {
            int client = parent2.chromT[i];
            if (!added[client]) {
                result.chromT.push_back(client);
                added[client] = true;
            }
        }
    }

    // Compléter l'individu avec l'algorithme Split
    split.generalSplit(result, parent1.eval.nbRoutes);
}

Genetic::Genetic(Params & params) :
        params(params),
        split(params),
        localSearch(params),
        population(params,this->split,this->localSearch),
        offspring(params){}
