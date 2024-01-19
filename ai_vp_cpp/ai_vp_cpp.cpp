// ai_vp_cpp.cpp: 定義應用程式的進入點。
//

#ifdef _MSC_VER
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#endif

#include "ai_vp_cpp.h"
#include "random_forest_200.h"

//
// Copyright (c) 2016-2019 Vinnie Falco (vinnie dot falco at gmail dot com)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Official repository: https://github.com/boostorg/beast
//

//------------------------------------------------------------------------------
//
// Example: Advanced server
//
//------------------------------------------------------------------------------

#include <boost/beast/core.hpp>
#include <boost/beast/http.hpp>
#include <boost/beast/websocket.hpp>
#include <boost/beast/version.hpp>
#include <boost/asio/bind_executor.hpp>
#include <boost/asio/dispatch.hpp>
#include <boost/asio/signal_set.hpp>
#include <boost/asio/strand.hpp>
#include <boost/make_unique.hpp>
#include <boost/optional.hpp>
#include <algorithm>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <thread>
#include <vector>
#include <unordered_map>
#include <charconv>
#include <regex>
#include <simdjson.h>
#include <nlohmann/json.hpp>
#include <malloc.h>

namespace beast = boost::beast;                 // from <boost/beast.hpp>
namespace http = beast::http;                   // from <boost/beast/http.hpp>
namespace websocket = beast::websocket;         // from <boost/beast/websocket.hpp>
namespace net = boost::asio;                    // from <boost/asio.hpp>
using tcp = boost::asio::ip::tcp;               // from <boost/asio/ip/tcp.hpp>

struct inputVariantInfo {
	std::string_view Chr;
	std::string_view Start;
	std::string_view Ref;
	std::string_view Alt;
	std::string_view Gene;
	std::string_view Otherinfo1;
	std::string_view Otherinfo2;
	std::string_view Otherinfo3;
	std::string_view MaxAlleleFreq;
	std::string_view SIFT;
	std::string_view PolyphenHDIV;
	std::string_view PolyphenHVAR;
	std::string_view LRT;
	std::string_view MutationTaster;
	std::string_view MutationAssessor;
	std::string_view FATHMM;
	std::string_view PROVEAN;
	std::string_view MetaSVM;
	std::string_view MetaLR;
	std::string_view MCAP;
	std::string_view CADD;
	std::string_view fathmmMKL;
	std::string_view Inheritance;
	std::string_view Consequence;
	std::string_view MaxEntScan;
	std::string_view FuncRefgene;
	std::string_view HGMD;
	std::string_view ClinVar;
};

struct variantInfo {
	float Rank;
	float Otherinfo1;
	float Otherinfo2;
	float Otherinfo3;
	float MaxAlleleFreq;
	float SIFT_converted_rankscore;
	float Polyphen2_HDIV_rankscore;
	float Polyphen2_HVAR_rankscore;
	float LRT_converted_rankscore;
	float MutationTaster_converted_rankscore;
	float MutationAssessor_score_rankscore;
	float FATHMM_converted_rankscore;
	float PROVEAN_converted_rankscore;
	float MetaSVM_rankscore;
	float MetaLR_rankscore;
	float M_CAP_rankscore;
	float CADD_raw_rankscore;
	float fathmm_MKL_coding_rankscore;
	float Consequence_3_prime_UTR_variant;
	float Consequence_5_prime_UTR_variant;
	float Consequence_coding_sequence_variant;
	float Consequence_downstream_gene_variant;
	float Consequence_frameshift_variant;
	float Consequence_inframe_deletion;
	float Consequence_inframe_insertion;
	float Consequence_intron_variant;
	float Consequence_missense_variant;
	float Consequence_non_coding_transcript_exon_variant;
	float Consequence_regulatory_region_variant;
	float Consequence_splice_acceptor_variant;
	float Consequence_splice_donor_variant;
	float Consequence_splice_region_variant;
	float Consequence_start_lost;
	float Consequence_stop_gained;
	float Consequence_stop_lost;
	float Consequence_stop_retained_variant;
	float Consequence_synonymous_variant;
	float Consequence_upstream_gene_variant;
	float MaxEntScan_Significance_yes;
	float Inheritance_AR;
	float Inheritance_S;
	float Inheritance_AD;
	float Inheritance_Multi;
	float Inheritance_QTL;
	float Inheritance_Mi;
	float Inheritance_IC;
	float Inheritance_XL;
	float Inheritance_XLR;
	float Inheritance_XLD;
	float Func_refgene_exonic;
	float Func_refgene_intronic;
	float Func_refgene_splicing;
	float HGMD_DM;
	float HGMD_DFP;
	float HGMD_DP;
	float HGMD_DM_;
	float HGMD_FP;
	float HGMD_R;
	float ClinVar_Conflicting;
	float ClinVar_UncertainSignificance;
	float ClinVar_LikelyBenign;
	float ClinVar_Benign;
	float ClinVar_Pathogenic;
	float ClinVar_LikelyPathogenic;
	float ClinVar_protective;
	float ClinVar_Riskfactor;
	float ClinVar_association;
	float ClinVar_drug_response;
};

struct resultInfo {
	uint64_t ModelRank;
	float Probability;
	uint64_t GenePrioritizer;
	std::string_view Gene;
	std::string_view Chr;
	std::string_view Start;
	std::string_view Ref;
	std::string_view Alt;
};

// For encoding/decoding to/from json, you need to provide a to_json method
static void to_json(nlohmann::json& j, const resultInfo& r) {
	j = nlohmann::json{ {"ModelRank", r.ModelRank}, {"Probability", r.Probability}, {"GenePrioritizer", r.GenePrioritizer}, {"Gene", r.Gene}, {"Chr", r.Chr}, {"Start", r.Start}, {"Ref", r.Ref}, {"Alt", r.Alt} };
}

static std::string performHttpPostRequest(
	const std::string& server,
	const std::string& port,
	const std::string& path,
	const std::string& requestBody
) {
	try {
		// Create an io_context
		boost::asio::io_context io_context;

		// Resolve the IP address and port
		tcp::resolver resolver(io_context);
		auto const results = resolver.resolve(server, port);

		// Create and connect the socket
		tcp::socket socket(io_context);
		boost::asio::connect(socket, results.begin(), results.end());

		// Set up the HTTP request
		http::request<http::string_body> req(http::verb::post, path, 11);
		req.set(http::field::host, server + ":" + port);
		req.set(http::field::content_type, "application/x-www-form-urlencoded; charset=UTF-8");

		// Build the request payload
		req.body() = requestBody;
		req.prepare_payload();

		// Send the HTTP request
		http::write(socket, req);

		// This buffer is used for reading and must be persisted
		beast::flat_buffer buffer;

		// Declare a container to hold the response
		http::response<http::dynamic_body> res;

		// Receive the HTTP response
		http::read(socket, buffer, res);

		// Check if the response is OK
		if (res.result() == http::status::ok || res.result() == http::status::created) {
			// Access the response body and convert it to a string
			std::string responseBody = boost::beast::buffers_to_string(res.body().data());

			// Parse the JSON response or process the string as needed
			return responseBody;
		}
		else {
			std::cerr << "HTTP response error: " << res.result_int() << std::endl;
			return std::string("");
		}
	}
	catch (const std::exception& e) {
		std::cerr << "Error: " << e.what() << std::endl;
		return std::string("");
	}
}

static auto METAMAP_HOST = std::getenv("METAMAP_HOST");
static auto METAMAP_PORT = std::getenv("METAMAP_PORT");
static auto METAMAP_ENDPOINT = std::getenv("METAMAP_ENDPOINT");

static std::string getKeywordFromMetamap(const std::string& abstract) {
	// Build the request payload
	auto payload = "abstract=" + abstract + "&tool=MetaMap";

	// Perform the HTTP request
	return performHttpPostRequest(METAMAP_HOST, METAMAP_PORT, METAMAP_ENDPOINT, payload);
}

static auto VP_HOST = std::getenv("VP_HOST");
static auto VP_PORT = std::getenv("VP_PORT");
static auto VP_ENDPOINT = std::getenv("VP_ENDPOINT");

static std::string getRankFromVariantPrioritizer(const std::string& keywords) {
	// Build the request payload
	auto payload = "genes=&keywords=" + keywords;

	// Perform the HTTP request
	return performHttpPostRequest(VP_HOST, VP_PORT, VP_ENDPOINT, payload);
}

// Report a failure
static void
fail(beast::error_code ec, char const* what)
{
	std::cerr << what << ": " << ec.message() << "\n";
}

// Handles an HTTP server connection
class http_session : public std::enable_shared_from_this<http_session>
{
	// This queue is used for HTTP pipelining.
	class queue
	{
		enum
		{
			// Maximum number of responses we will queue
			limit = 8
		};

		// The type-erased, saved work item
		struct work
		{
			virtual ~work() = default;
			virtual void operator()() = 0;
		};

		http_session& self_;
		std::vector<std::unique_ptr<work>> items_;

	public:
		explicit
			queue(http_session& self)
			: self_(self)
		{
			static_assert(limit > 0, "queue limit must be positive");
			items_.reserve(limit);
		}

		// Returns `true` if we have reached the queue limit
		bool
			is_full() const
		{
			return items_.size() >= limit;
		}

		// Called when a message finishes sending
		// Returns `true` if the caller should initiate a read
		bool
			on_write()
		{
			BOOST_ASSERT(!items_.empty());
			auto const was_full = is_full();
			items_.erase(items_.begin());
			if (!items_.empty())
				(*items_.front())();
			return was_full;
		}

		// Called by the HTTP handler to send a response.
		template<bool isRequest, class Body, class Fields>
		void
			operator()(http::message<isRequest, Body, Fields>&& msg)
		{
			// This holds a work item
			struct work_impl : work
			{
				http_session& self_;
				http::message<isRequest, Body, Fields> msg_;

				work_impl(
					http_session& self,
					http::message<isRequest, Body, Fields>&& msg)
					: self_(self)
					, msg_(std::move(msg))
				{
				}

				void
					operator()()
				{
					http::async_write(
						self_.stream_,
						msg_,
						beast::bind_front_handler(
							&http_session::on_write,
							self_.shared_from_this(),
							msg_.need_eof()));
				}
			};

			// Allocate and store the work
			items_.push_back(
				boost::make_unique<work_impl>(self_, std::move(msg)));

			// If there was no previous work, start this one
			if (items_.size() == 1)
				(*items_.front())();
		}
	};

	beast::tcp_stream stream_;
	beast::flat_buffer buffer_;
	queue queue_;

	// The parser is stored in an optional container so we can
	// construct it from scratch it at the beginning of each new message.
	boost::optional<http::request_parser<http::buffer_body>> parser_;

public:
	// Take ownership of the socket
	http_session(
		tcp::socket&& socket)
		: stream_(std::move(socket))
		, queue_(*this)
	{
	}

	// Start the session
	void
		run()
	{
		// We need to be executing within a strand to perform async operations
		// on the I/O objects in this session. Although not strictly necessary
		// for single-threaded contexts, this example code is written to be
		// thread-safe by default.
		net::dispatch(
			stream_.get_executor(),
			beast::bind_front_handler(
				&http_session::do_read,
				this->shared_from_this()));
	}

private:
	void
		do_read()
	{
		// Construct a new parser for each message
		parser_.emplace();

		// Apply a reasonable limit to the allowed size
		// of the body in bytes to prevent abuse.
		parser_->body_limit(std::numeric_limits<uint64_t>().max());

		// Set the timeout.
		//stream_.expires_after(std::chrono::seconds(30));

		// Read a request using the parser-oriented interface
		http::async_read(
			stream_,
			buffer_,
			*parser_,
			beast::bind_front_handler(
				&http_session::on_read,
				shared_from_this()));
	}

	void
		on_read(beast::error_code ec, std::size_t bytes_transferred)
	{
#ifdef __GLIBC__
		malloc_trim(0);
#endif

		boost::ignore_unused(bytes_transferred);

		// This means they closed the connection
		if (ec == http::error::end_of_stream)
			return do_close();

		auto& p = parser_.value();

		http::read_header(stream_, buffer_, p, ec);
		if (ec) {
			return fail(ec, "read header");
		}

		try {
			std::stringstream ss;
			while (!p.is_done())
			{
				char buf[8192]{};
				p.get().body().data = buf;
				p.get().body().size = sizeof(buf);
				http::read(stream_, buffer_, p, ec);
				if (ec == http::error::need_buffer) {
					ec = {};
				}

				if (ec) {
					return fail(ec, "read body");
				}

				ss.write(buf, sizeof(buf) - p.get().body().size);
			}

			auto padding_json = simdjson::padded_string(ss.str());
			simdjson::ondemand::parser parser;
			auto doc = parser.iterate(padding_json);
			auto results = std::vector<resultInfo>();
			bool isMetadata = true;
			auto gene_rank_dict = std::unordered_map<std::string, uint64_t>();
			for (auto e : doc) {
				if (isMetadata) {
					auto keywords = std::string(e["keywords"].get_string().value());
					auto isAbstract = e["isAbstract"].get_bool().value();

					if (isAbstract) {
						auto res = getKeywordFromMetamap(keywords);

						try {
							auto padding_json = simdjson::padded_string(res);
							simdjson::ondemand::parser parser;
							auto doc = parser.iterate(padding_json);

							keywords = std::string(doc["keyword"].get_string().value());
						}
						catch (const std::exception& e) {
							std::cerr << "Error parsing JSON: " << e.what() << std::endl;
							keywords = "";
						}
					}

					auto res = getRankFromVariantPrioritizer(keywords);
					if (res != std::string(""))
					{
						try {
							auto padding_json = simdjson::padded_string(res);
							simdjson::ondemand::parser parser;
							auto doc = parser.iterate(padding_json);

							for (auto field : doc) {
								auto gene = std::string(field["geneName"].get_string().value());
								auto rank = field["rank"].get_uint64().value();
								gene_rank_dict[gene] = rank;
							}
						}
						catch (const std::exception& e) {
							std::cerr << "Error parsing JSON: " << e.what() << std::endl;
						}
					}

					isMetadata = false;
				}
				else {
					auto i = inputVariantInfo();
					i.Chr = e["Chr"].is_null() ? "" : e["Chr"].get_string().value();
					i.Start = e["Start"].is_null() ? "" : e["Start"].get_string().value();
					i.Ref = e["Ref"].is_null() ? "" : e["Ref"].get_string().value();
					i.Alt = e["Alt"].is_null() ? "" : e["Alt"].get_string().value();
					i.Gene = e["Gene"].is_null() ? "" : e["Gene"].get_string().value();
					i.Otherinfo1 = e["Otherinfo1"].is_null() ? "" : e["Otherinfo1"].get_string().value();
					i.Otherinfo2 = e["Otherinfo2"].is_null() ? "" : e["Otherinfo2"].get_string().value();
					i.Otherinfo3 = e["Otherinfo3"].is_null() ? "" : e["Otherinfo3"].get_string().value();
					i.MaxAlleleFreq = e["MaxAlleleFreq"].is_null() ? "" : e["MaxAlleleFreq"].get_string().value();
					i.SIFT = e["SIFT"].is_null() ? "" : e["SIFT"].get_string().value();
					i.PolyphenHDIV = e["PolyphenHDIV"].is_null() ? "" : e["PolyphenHDIV"].get_string().value();
					i.PolyphenHVAR = e["PolyphenHVAR"].is_null() ? "" : e["PolyphenHVAR"].get_string().value();
					i.LRT = e["LRT"].is_null() ? "" : e["LRT"].get_string().value();
					i.MutationTaster = e["MutationTaster"].is_null() ? "" : e["MutationTaster"].get_string().value();
					i.MutationAssessor = e["MutationAssessor"].is_null() ? "" : e["MutationAssessor"].get_string().value();
					i.FATHMM = e["FATHMM"].is_null() ? "" : e["FATHMM"].get_string().value();
					i.PROVEAN = e["PROVEAN"].is_null() ? "" : e["PROVEAN"].get_string().value();
					i.MetaSVM = e["MetaSVM"].is_null() ? "" : e["MetaSVM"].get_string().value();
					i.MetaLR = e["MetaLR"].is_null() ? "" : e["MetaLR"].get_string().value();
					i.MCAP = e["MCAP"].is_null() ? "" : e["MCAP"].get_string().value();
					i.CADD = e["CADD"].is_null() ? "" : e["CADD"].get_string().value();
					i.fathmmMKL = e["fathmmMKL"].is_null() ? "" : e["fathmmMKL"].get_string().value();
					i.Inheritance = e["Inheritance"].is_null() ? "" : e["Inheritance"].get_string().value();
					i.Consequence = e["Consequence"].is_null() ? "" : e["Consequence"].get_string().value();
					i.MaxEntScan = e["MaxEntScan"].is_null() ? "" : e["MaxEntScan"].get_string().value();
					i.FuncRefgene = e["FuncRefgene"].is_null() ? "" : e["FuncRefgene"].get_string().value();
					i.HGMD = e["HGMD"].is_null() ? "" : e["HGMD"].get_string().value();
					i.ClinVar = e["ClinVar"].is_null() ? "" : e["ClinVar"].get_string().value();

					auto v = variantInfo{};

					auto delimiter = ";";
					auto regexDelimiter = std::regex(delimiter);
					auto text = std::string(i.Gene);
					auto it = std::sregex_token_iterator(text.begin(), text.end(), regexDelimiter, -1);
					auto end = std::sregex_token_iterator();

					auto genes = std::vector<std::string>();
					while (it != end) {
						std::string token = it->str();

						genes.emplace_back(token);

						++it;
					}

					auto max_gene_ptr = std::max_element(genes.begin(), genes.end(), [&](const auto& a, const auto& b) {
						return gene_rank_dict[a] < gene_rank_dict[b];
						});

					auto rank = gene_rank_dict[*max_gene_ptr];
					v.Rank = static_cast<float>(rank);

					if (i.Otherinfo1 == "het") { v.Otherinfo1 = 1; }

					auto to_float = [](std::string_view s) {
						float value{};
#if __cpp_lib_to_chars >= 202306L
						if (std::from_chars(s.data(), s.data() + s.size(), value))
#else
						if (std::from_chars(s.data(), s.data() + s.size(), value).ec == std::errc{})
#endif
							return value;
						else
							return std::numeric_limits<float>::quiet_NaN();
						};

					v.Otherinfo2 = to_float(i.Otherinfo2);
					if (std::isnan(v.Otherinfo2)) { v.Otherinfo2 = 0; }

					v.Otherinfo3 = to_float(i.Otherinfo3);
					if (std::isnan(v.Otherinfo3)) { v.Otherinfo3 = 0; }

					v.MaxAlleleFreq = to_float(i.MaxAlleleFreq);
					if (std::isnan(v.MaxAlleleFreq)) { v.MaxAlleleFreq = 0; }

					v.SIFT_converted_rankscore = to_float(i.SIFT);
					if (std::isnan(v.SIFT_converted_rankscore)) { v.SIFT_converted_rankscore = -1; }

					v.Polyphen2_HDIV_rankscore = to_float(i.PolyphenHDIV);
					if (std::isnan(v.Polyphen2_HDIV_rankscore)) { v.Polyphen2_HDIV_rankscore = -1; }

					v.Polyphen2_HVAR_rankscore = to_float(i.PolyphenHVAR);
					if (std::isnan(v.Polyphen2_HVAR_rankscore)) { v.Polyphen2_HVAR_rankscore = -1; }

					v.LRT_converted_rankscore = to_float(i.LRT);
					if (std::isnan(v.LRT_converted_rankscore)) { v.LRT_converted_rankscore = -1; }

					v.MutationTaster_converted_rankscore = to_float(i.MutationTaster);
					if (std::isnan(v.MutationTaster_converted_rankscore)) { v.MutationTaster_converted_rankscore = -1; }

					v.MutationAssessor_score_rankscore = to_float(i.MutationAssessor);
					if (std::isnan(v.MutationAssessor_score_rankscore)) { v.MutationAssessor_score_rankscore = -1; }

					v.FATHMM_converted_rankscore = to_float(i.FATHMM);
					if (std::isnan(v.FATHMM_converted_rankscore)) { v.FATHMM_converted_rankscore = -1; }

					v.PROVEAN_converted_rankscore = to_float(i.PROVEAN);
					if (std::isnan(v.PROVEAN_converted_rankscore)) { v.PROVEAN_converted_rankscore = -1; }

					v.MetaSVM_rankscore = to_float(i.MetaSVM);
					if (std::isnan(v.MetaSVM_rankscore)) { v.MetaSVM_rankscore = -1; }

					v.MetaLR_rankscore = to_float(i.MetaLR);
					if (std::isnan(v.MetaLR_rankscore)) { v.MetaLR_rankscore = -1; }

					v.M_CAP_rankscore = to_float(i.MCAP);
					if (std::isnan(v.M_CAP_rankscore)) { v.M_CAP_rankscore = -1; }

					v.CADD_raw_rankscore = to_float(i.CADD);
					if (std::isnan(v.CADD_raw_rankscore)) { v.CADD_raw_rankscore = -1; }

					v.fathmm_MKL_coding_rankscore = to_float(i.fathmmMKL);
					if (std::isnan(v.fathmm_MKL_coding_rankscore)) { v.fathmm_MKL_coding_rankscore = -1; }

					delimiter = ",";
					regexDelimiter = std::regex(delimiter);
					text = std::string(i.Consequence);
					it = std::sregex_token_iterator(text.begin(), text.end(), regexDelimiter, -1);
					end = std::sregex_token_iterator();

					while (it != end) {
						std::string token = it->str();

						if (token == "3_prime_UTR_variant") { v.Consequence_3_prime_UTR_variant = 1; }
						else if (token == "5_prime_UTR_variant") { v.Consequence_5_prime_UTR_variant = 1; }
						else if (token == "coding_sequence_variant") { v.Consequence_coding_sequence_variant = 1; }
						else if (token == "downstream_gene_variant") { v.Consequence_downstream_gene_variant = 1; }
						else if (token == "frameshift_variant") { v.Consequence_frameshift_variant = 1; }
						else if (token == "inframe_deletion") { v.Consequence_inframe_deletion = 1; }
						else if (token == "inframe_insertion") { v.Consequence_inframe_insertion = 1; }
						else if (token == "intron_variant") { v.Consequence_intron_variant = 1; }
						else if (token == "missense_variant") { v.Consequence_missense_variant = 1; }
						else if (token == "non_coding_transcript_exon_variant") { v.Consequence_non_coding_transcript_exon_variant = 1; }
						else if (token == "regulatory_region_variant") { v.Consequence_regulatory_region_variant = 1; }
						else if (token == "splice_acceptor_variant") { v.Consequence_splice_acceptor_variant = 1; }
						else if (token == "splice_donor_variant") { v.Consequence_splice_donor_variant = 1; }
						else if (token == "splice_region_variant") { v.Consequence_splice_region_variant = 1; }
						else if (token == "start_lost") { v.Consequence_start_lost = 1; }
						else if (token == "stop_gained") { v.Consequence_stop_gained = 1; }
						else if (token == "stop_lost") { v.Consequence_stop_lost = 1; }
						else if (token == "stop_retained_variant") { v.Consequence_stop_retained_variant = 1; }
						else if (token == "synonymous_variant") { v.Consequence_synonymous_variant = 1; }
						else if (token == "upstream_gene_variant") { v.Consequence_upstream_gene_variant = 1; }

						++it;
					}

					if (i.MaxEntScan == "yes") { v.MaxEntScan_Significance_yes = 1; }

					delimiter = ", ";
					regexDelimiter = std::regex(delimiter);
					text = std::string(i.Inheritance);
					it = std::sregex_token_iterator(text.begin(), text.end(), regexDelimiter, -1);
					end = std::sregex_token_iterator();

					while (it != end) {
						std::string token = it->str();

						if (token == "AR") { v.Inheritance_AR = 1; }
						else if (token == "S") { v.Inheritance_S = 1; }
						else if (token == "AD") { v.Inheritance_AD = 1; }
						else if (token == "Multi") { v.Inheritance_Multi = 1; }
						else if (token == "QTL") { v.Inheritance_QTL = 1; }
						else if (token == "Mi") { v.Inheritance_Mi = 1; }
						else if (token == "IC") { v.Inheritance_IC = 1; }
						else if (token == "XL") { v.Inheritance_XL = 1; }
						else if (token == "XLR") { v.Inheritance_XLR = 1; }
						else if (token == "XLD") { v.Inheritance_XLD = 1; }

						++it;
					}

					delimiter = ";";
					regexDelimiter = std::regex(delimiter);
					text = std::string(i.FuncRefgene);
					it = std::sregex_token_iterator(text.begin(), text.end(), regexDelimiter, -1);
					end = std::sregex_token_iterator();

					while (it != end) {
						std::string token = it->str();

						if (token == "exonic") { v.Func_refgene_exonic = 1; }
						else if (token == "intronic") { v.Func_refgene_intronic = 1; }
						else if (token == "splicing") { v.Func_refgene_splicing = 1; }

						++it;
					}

					delimiter = ";";
					regexDelimiter = std::regex(delimiter);
					text = std::string(i.HGMD);
					it = std::sregex_token_iterator(text.begin(), text.end(), regexDelimiter, -1);
					end = std::sregex_token_iterator();

					while (it != end) {
						std::string token = it->str();

						if (token == "DM") { v.HGMD_DM = 1; }
						else if (token == "DFP") { v.HGMD_DFP = 1; }
						else if (token == "DP") { v.HGMD_DP = 1; }
						else if (token == "DM?") { v.HGMD_DM_ = 1; }
						else if (token == "FP") { v.HGMD_FP = 1; }
						else if (token == "R") { v.HGMD_R = 1; }

						++it;
					}

					delimiter = "; ";
					regexDelimiter = std::regex(delimiter);
					text = std::string(i.ClinVar);
					it = std::sregex_token_iterator(text.begin(), text.end(), regexDelimiter, -1);
					end = std::sregex_token_iterator();

					while (it != end) {
						std::string token = it->str();

						if (token == "Conflicting interpretations of pathogenicity") { v.ClinVar_Conflicting = 1; }
						else if (token == "Uncertain significance") { v.ClinVar_UncertainSignificance = 1; }
						else if (token == "Likely benign") { v.ClinVar_LikelyBenign = 1; }
						else if (token == "Benign") { v.ClinVar_Benign = 1; }
						else if (token == "Pathogenic") { v.ClinVar_Pathogenic = 1; }
						else if (token == "Likely pathogenic") { v.ClinVar_LikelyPathogenic = 1; }
						else if (token == "protective") { v.ClinVar_protective = 1; }
						else if (token == "risk factor") { v.ClinVar_Riskfactor = 1; }
						else if (token == "association") { v.ClinVar_association = 1; }
						else if (token == "drug response") { v.ClinVar_drug_response = 1; }

						++it;
					}

					auto floatArray = reinterpret_cast<float*>(&v);
					auto result = predict(floatArray);
					auto floatValue = static_cast<float>(result) / 100.0f;

					auto r = resultInfo();
					r.ModelRank = 0;
					r.Probability = floatValue;
					r.GenePrioritizer = rank;
					r.Gene = i.Gene;
					r.Chr = i.Chr;
					r.Start = i.Start;
					r.Ref = i.Ref;
					r.Alt = i.Alt;

					results.emplace_back(r);
				}
			}

			std::sort(results.begin(), results.end(), [](const auto& a, const auto& b) {
				return a.Probability > b.Probability;
				});

			for (size_t i = 0; i < results.size(); ++i) {
				results[i].ModelRank = i + 1;
			}

			nlohmann::json json_array = results;

			http::response<http::string_body> res{ http::status::ok, 11 };
			res.set(http::field::content_type, "application/json");
			res.body() = json_array.dump();
			res.prepare_payload();
			queue_(std::move(res));
		}
		catch (const std::exception& e) {
			std::cerr << "Error: " << e.what() << std::endl;
			http::response<http::string_body> res{ http::status::internal_server_error, 11 };
			res.set(http::field::content_type, "text/plain");
			res.body() = "An error occurred.";
			res.prepare_payload();
			queue_(std::move(res));
		}

#ifdef __GLIBC__
		malloc_trim(0);
#endif

		// If we aren't at the queue limit, try to pipeline another request
		if (!queue_.is_full())
			do_read();
	}

	void
		on_write(bool close, beast::error_code ec, std::size_t bytes_transferred)
	{
		boost::ignore_unused(bytes_transferred);

		if (ec)
			return fail(ec, "write");

		if (close)
		{
			// This means we should close the connection, usually because
			// the response indicated the "Connection: close" semantic.
			return do_close();
		}

		// Inform the queue that a write completed
		if (queue_.on_write())
		{
			// Read another request
			do_read();
		}
	}

	void
		do_close()
	{
		// Send a TCP shutdown
		beast::error_code ec;
		stream_.socket().shutdown(tcp::socket::shutdown_send, ec);

		// At this point the connection is closed gracefully
	}
};

//------------------------------------------------------------------------------

// Accepts incoming connections and launches the sessions
class listener : public std::enable_shared_from_this<listener>
{
	net::io_context& ioc_;
	tcp::acceptor acceptor_;

public:
	listener(
		net::io_context& ioc,
		tcp::endpoint endpoint)
		: ioc_(ioc)
		, acceptor_(net::make_strand(ioc))
	{
		beast::error_code ec;

		// Open the acceptor
		acceptor_.open(endpoint.protocol(), ec);
		if (ec)
		{
			fail(ec, "open");
			return;
		}

		// Allow address reuse
		acceptor_.set_option(net::socket_base::reuse_address(true), ec);
		if (ec)
		{
			fail(ec, "set_option");
			return;
		}

		// Bind to the server address
		acceptor_.bind(endpoint, ec);
		if (ec)
		{
			fail(ec, "bind");
			return;
		}

		// Start listening for connections
		acceptor_.listen(
			net::socket_base::max_listen_connections, ec);
		if (ec)
		{
			fail(ec, "listen");
			return;
		}
	}

	// Start accepting incoming connections
	void
		run()
	{
		// We need to be executing within a strand to perform async operations
		// on the I/O objects in this session. Although not strictly necessary
		// for single-threaded contexts, this example code is written to be
		// thread-safe by default.
		net::dispatch(
			acceptor_.get_executor(),
			beast::bind_front_handler(
				&listener::do_accept,
				this->shared_from_this()));
	}

private:
	void
		do_accept()
	{
		// The new connection gets its own strand
		acceptor_.async_accept(
			net::make_strand(ioc_),
			beast::bind_front_handler(
				&listener::on_accept,
				shared_from_this()));
	}

	void
		on_accept(beast::error_code ec, tcp::socket socket)
	{
		if (ec)
		{
			fail(ec, "accept");
		}
		else
		{
			// Create the http session and run it
			std::make_shared<http_session>(
				std::move(socket))->run();
		}

		// Accept another connection
		do_accept();
	}
};

//------------------------------------------------------------------------------

int main(int argc, char* argv[])
{
	auto const address = net::ip::make_address("0.0.0.0");
	auto const port = static_cast<unsigned short>(8080);
	auto const threads = std::max<int>(1, std::thread::hardware_concurrency());

	// The io_context is required for all I/O
	net::io_context ioc{ threads };

	// Create and launch a listening port
	std::make_shared<listener>(
		ioc,
		tcp::endpoint{ address, port })->run();

	// Capture SIGINT and SIGTERM to perform a clean shutdown
	net::signal_set signals(ioc, SIGINT, SIGTERM);
	signals.async_wait(
		[&](beast::error_code const&, int)
		{
			// Stop the `io_context`. This will cause `run()`
			// to return immediately, eventually destroying the
			// `io_context` and all of the sockets in it.
			ioc.stop();
		});

	// Run the I/O service on the requested number of threads
	std::vector<std::thread> v;
	v.reserve(threads - 1);
	for (auto i = threads - 1; i > 0; --i)
		v.emplace_back(
			[&ioc]
			{
				ioc.run();
			});
	ioc.run();

	// (If we get here, it means we got a SIGINT or SIGTERM)

	// Block until all the threads exit
	for (auto& t : v)
		t.join();

#ifdef _MSC_VER
	_CrtDumpMemoryLeaks();
#endif

	return EXIT_SUCCESS;
}
