using System;
using System.Linq;
using System.Threading.Tasks;

class GeneticAlgorithm
{
    static Random rand = new Random();
    const int PopulationSize = 100;
    const int GeneCount = 1; // Assuming a single gene for simplicity
    const int MaxGenerations = 100;
    const double MutationRate = 0.01;
    const double CrossoverRate = 0.7;
    const double MinValue = -10.0;
    const double MaxValue = 10.0;

    static void Main(string[] args)
    {
        // Initialize population
        double[][] population = new double[PopulationSize][];
        for (int i = 0; i < PopulationSize; i++)
        {
            population[i] = new double[GeneCount];
            for (int j = 0; j < GeneCount; j++)
            {
                population[i][j] = rand.NextDouble() * (MaxValue - MinValue) + MinValue;
            }
        }

        int generation = 0;
        while (generation < MaxGenerations)
        {
            // Evaluate fitness in parallel
            double[] fitness = new double[PopulationSize];
            Parallel.For(0, PopulationSize, i =>
            {
                fitness[i] = FitnessFunction(population[i][0]);
            });

            // Selection
            double[][] selected = TournamentSelection(population, fitness);

            // Crossover
            double[][] offspring = Crossover(selected);

            // Mutation
            Mutate(offspring);

            // Replacement with Elitism
            int bestIndex = Array.IndexOf(fitness, fitness.Max());
            population = offspring;
            population[0] = population[bestIndex]; // Elitism: carry the best to the next generation

            // Print best solution every generation
            Console.WriteLine($"Generation {generation}, Best Fitness: {fitness[bestIndex]:F2}, x: {population[bestIndex][0]:F2}");

            generation++;
        }
    }

    static double FitnessFunction(double x) => x * x;

    static double[][] TournamentSelection(double[][] population, double[] fitness)
    {
        double[][] selected = new double[PopulationSize][];
        for (int i = 0; i < PopulationSize; i++)
        {
            // Randomly select two individuals from the population
            int first = rand.Next(PopulationSize);
            int second = rand.Next(PopulationSize);
            // Choose the one with better fitness to be part of the mating pool
            selected[i] = fitness[first] > fitness[second] ? population[first] : population[second];
        }
        return selected;
    }


    static double[][] Crossover(double[][] selected)
    {
        double[][] offspring = new double[PopulationSize][];
        for (int i = 0; i < PopulationSize; i += 2)
        {
            offspring[i] = new double[GeneCount];
            offspring[i + 1] = new double[GeneCount];
            for (int j = 0; j < GeneCount; j++)
            {
                // Perform crossover with a certain probability
                if (rand.NextDouble() < CrossoverRate)
                {
                    // If crossover occurs, offspring inherit genes directly from parents
                    offspring[i][j] = selected[i][j];
                    offspring[i + 1][j] = selected[i + 1][j];
                }
                else
                {
                    // If crossover does not occur, swap genes between parents
                    offspring[i][j] = selected[i + 1][j];
                    offspring[i + 1][j] = selected[i][j];
                }
            }
        }
        return offspring;
    }

    static void Mutate(double[][] offspring)
    {
        for (int i = 0; i < offspring.Length; i++)
        {
            for (int j = 0; j < GeneCount; j++)
            {
                // With a certain probability, mutate a gene
                if (rand.NextDouble() < MutationRate)
                {
                    // Apply a small change to the gene value
                    offspring[i][j] += (rand.NextDouble() * 2 - 1) * (MaxValue - MinValue) * MutationRate;
                    // Ensure the mutated gene is within the predefined bounds
                    offspring[i][j] = Math.Max(MinValue, Math.Min(MaxValue, offspring[i][j]));
                }
            }
        }
    }
}
