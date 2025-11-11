static StatsOut computeStats(const std::vector<std::vector<std::vector<int>>>& statM) {
    const size_t I = statM.size(); const size_t M = I? statM[0].size() : 0; const size_t R = (M? statM[0][0].size() : 0);
    StatsOut out; out.mean.assign(M, std::vector<double>(R, 0.0)); out.relErr.assign(M, std::vector<double>(R, 0.0)); out.sum.assign(M, std::vector<int>(R, 0));
    for (size_t m=0; m<M; ++m) for (size_t r=0; r<R; ++r) {
        long long S=0; long double Q=0.0L;
        for (size_t i=0;i<I;++i){ int c=statM[i][m][r]; S+=c; Q+=1.0L*c*c; }
        const double N=double(I); const double mu=(I? double(S)/N : 0.0);
        const double var=(I>1)? double((Q - (1.0L*S*S)/N)/(N-1.0)) : 0.0;
        const double se=(I>0)? std::sqrt(std::max(0.0, var)/N) : 0.0;
        out.sum[m][r]=int(S); out.mean[m][r]=mu; out.relErr[m][r]=(mu>0.0)? se/mu : 0.0;
    }
    return out;
}
static void printStatsOut(const StatsOut& S, const std::vector<std::string>& matNames, const std::vector<int>& MTs, std::ostream& os = std::cout) {
    const size_t M=S.sum.size(); if (!M) { os << "(no data)\n"; return; } const size_t R=S.sum[0].size();
    std::vector<long long> rowTot(M,0), colTot(R,0); long long grand=0;
    for (size_t i=0;i<M;++i) for (size_t j=0;j<R;++j){ long long v=S.sum[i][j]; rowTot[i]+=v; colTot[j]+=v; grand+=v; }
    auto pct=[&](long long x){ return grand? 100.0*double(x)/double(grand) : 0.0; };
    size_t wName=4; for (size_t i=0;i<std::min(M,matNames.size());++i) wName = std::max(wName, matNames[i].size());
    std::vector<size_t> wCol(R, 14);
    for (size_t j=0;j<R;++j){ wCol[j] = std::max(wCol[j], mtLabel(MTs[j]).size()); for (size_t i=0;i<M;++i){
        if (S.sum[i][j]==0) continue; std::ostringstream ss; double rPct = std::isfinite(S.relErr[i][j]) ? 100.0*S.relErr[i][j] : 0.0;
        ss << std::setprecision(3) << std::scientific << S.mean[i][j] << " ± " << std::fixed << std::setprecision(1) << rPct << "% (" << S.sum[i][j] << ")"; wCol[j] = std::max<size_t>(wCol[j], ss.str().size());
    } }
    os << "\n=== Statistics ===\nTotal events: " << grand << "\n\n-- By reaction (counts) --\n";
    for (size_t j=0;j<R;++j) os << std::left << std::setw(int(wCol[j])) << mtLabel(MTs[j]) << "  " << std::right << std::setw(10) << colTot[j] << "  " << std::fixed << std::setprecision(2) << std::setw(6) << pct(colTot[j]) << "%\n";
    os << "\n-- By material (counts) --\n";
    for (size_t i=0;i<M;++i){ const std::string& name = (i<matNames.size()? matNames[i] : ("mat"+std::to_string(i)));
        os << std::left << std::setw(int(wName)) << name << "  " << std::right << std::setw(10) << rowTot[i] << "  " << std::fixed << std::setprecision(2) << std::setw(6) << pct(rowTot[i]) << "%\n"; }
    os << "\n-- Matrix: mean +- relErr% (count) --\n" << std::left << std::setw(int(wName)) << "" << "  ";
    for (size_t j=0;j<R;++j) os << std::left << std::setw(int(wCol[j])) << mtLabel(MTs[j]) << "  "; os << "\n";
    for (size_t i=0;i<M;++i){ const std::string& name = (i<matNames.size()? matNames[i] : ("mat"+std::to_string(i)));
        os << std::left << std::setw(int(wName)) << name << "  ";
        for (size_t j=0;j<R;++j){ if (S.sum[i][j]==0) os << std::left << std::setw(int(wCol[j])) << "-" << "  ";
            else { double rPct = std::isfinite(S.relErr[i][j]) ? 100.0*S.relErr[i][j] : 0.0; std::ostringstream cell;
                cell << std::setprecision(3) << std::scientific << S.mean[i][j] << " +- " << std::fixed << std::setprecision(1) << rPct << "% (" << S.sum[i][j] << ")";
                os << std::left << std::setw(int(wCol[j])) << cell.str() << "  "; } }
        os << "\n"; }
}