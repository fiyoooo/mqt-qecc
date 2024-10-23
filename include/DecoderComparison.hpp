#include "Code.hpp"

#include <tuple>
#include <vector>

class DecoderComparisonHelper {
public:
    DecoderComparisonHelper(Code& code, bool useMaxSATDecoder = false, bool usePeelDecoder = false)
        : code(code), pcm(code.toBpSparse()), maxsat(useMaxSATDecoder), peel(usePeelDecoder) {}

    /**
     * Tests all decoders with a given error.
     * @returns error estimation of each decoder and their respective runtimes
     */
    std::tuple<std::vector<gf2Vec>, std::vector<double>> testWithError(gf2Vec const& error, bool verbose = true);

    /**
     * Run multiple rounds of tests for all decoders with a given error probability.
     * @param prob probability of a bit-flip error
     * @param numRounds number of test rounds
     * @return error rates and average runtimes of each decoder
     */
    std::tuple<std::vector<double>, std::vector<double>> testWithErrorRate(double prob, int numRounds, bool verbose = true);

private:
    Code&              code;
    ldpc::bp::BpSparse pcm;
    bool               maxsat;
    bool               peel;

    static gf2Vec generateRandomBitFlipError(size_t size, double prob);

    std::tuple<std::vector<gf2Vec>, std::vector<double>> runDecoders(gf2Vec const& syndr);

    gf2Vec               runUFDecoder(gf2Vec const& syndrome);
    gf2Vec               runMaxSatDecoder(gf2Vec const& syndrome);
    std::vector<uint8_t> runUfdPeelDecoder(std::vector<uint8_t> const& syndrome);
    std::vector<uint8_t> runUfdMatrixDecoder(std::vector<uint8_t> const& syndrome);
    std::vector<uint8_t> runUfdMaxSatDecoder(std::vector<uint8_t> const& syndrome);

    template <typename Func, typename ReturnType = std::invoke_result_t<Func>>
    std::tuple<double, ReturnType> measureRuntime(Func func);

    bool isCorrectable(gf2Vec const& error, gf2Vec const& estimation);
};