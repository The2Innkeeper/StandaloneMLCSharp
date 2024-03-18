using System;
using System.Linq;

class GeneticAlgorithm
{
    static Random rand = new Random();
    const int PopulationSize = 100;
    const int GenomeSize = 32; // Increased for better precision
    const int MaxGenerations = 100;
    const double MutationRate = 0.01;
    const double CrossoverRate = 0.7;

    static void Main(string[] args)
    {
        // Initialize population
        bool[][] population = new bool[PopulationSize][];
        for (int i = 0; i < PopulationSize; i++)
        {
            population[i] = Enumerable.Range(0, GenomeSize).Select(_ => rand.NextDouble() >= 0.5).ToArray();
        }

        int generation = 0;
        while (generation < MaxGenerations)
        {
            // Evaluate fitness
            double[] fitness = population.Select(Decode).Select(FitnessFunction).ToArray();

            // Selection
            bool[][] selected = TournamentSelection(population, fitness);

            // Crossover
            bool[][] offspring = Crossover(selected);

            // Mutation
            Mutate(offspring);

            // Replacement
            population = offspring;

            // Print best solution every generation
            int maxFitnessIndex = Array.IndexOf(fitness, fitness.Max());
            Console.WriteLine($"Generation {generation}, Best Fitness: {fitness[maxFitnessIndex]:F2}, x: {Decode(population[maxFitnessIndex]):F2}");

            generation++;
        }
    }

    static double FitnessFunction(double x) => x * x;

    static double Decode(bool[] genome)
    {
        // Convert binary to decimal
        long value = 0;
        for (int i = 0; i < genome.Length; i++)
        {
            if (genome[i])
            {
                value |= 1L << i;
            }
        }

        // Map to range [-10, 10]
        double proportion = (double)value / (1L << GenomeSize);
        return proportion * 20.0 - 10.0;
    }

    static bool[][] TournamentSelection(bool[][] population, double[] fitness)
    {
        bool[][] selected = new bool[population.Length][];
        for (int i = 0; i < population.Length; i++)
        {
            int first = rand.Next(PopulationSize);
            int second = rand.Next(PopulationSize);
            selected[i] = fitness[first] > fitness[second] ? population[first] : population[second];
        }
        return selected;
    }

    static bool[][] Crossover(bool[][] selected)
    {
        bool[][] offspring = new bool[selected.Length][];
        for (int i = 0; i < selected.Length; i += 2)
        {
            offspring[i] = new bool[GenomeSize];
            offspring[i + 1] = new bool[GenomeSize];
            int crossoverPoint = rand.Next(GenomeSize);
            for (int j = 0; j < GenomeSize; j++)
            {
                if (j < crossoverPoint)
                {
                    offspring[i][j] = selected[i][j];
                    offspring[i + 1][j] = selected[i + 1][j];
                }
                else
                {
                    offspring[i][j] = selected[i + 1][j];
                    offspring[i + 1][j] = selected[i][j];
                }
            }
        }
        return offspring;
    }

    static void Mutate(bool[][] offspring)
    {
        for (int i = 0; i < offspring.Length; i++)
        {
            for (int j = 0; j < GenomeSize; j++)
            {
                if (rand.NextDouble() < MutationRate)
                {
                    offspring[i][j] = !offspring[i][j];
                }
            }
        }
    }
}
