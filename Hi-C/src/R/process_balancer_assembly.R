require(dplyr)
require(tidyr)

dir <- "."
assembly <- "3"
margin <- 1000L

as <- tbl_df(data.table::fread(paste0(dir, "/analysis/assembly_dm6bal", assembly, ".tab"), header = T, data.table = F)) %>%
  mutate(strand = ifelse(breakpoint2 > breakpoint1, "+", "-"))


as.track <- as %>%
  mutate(
    start = pmin(breakpoint1, breakpoint2) - 1L,
    end = pmax(breakpoint1, breakpoint2),
    score = 1000L,
    thickStart = start,
    thickEnd = end,
    itemRgb = apply(col2rgb(as$col), 2, paste, collapse = ",")
  ) %>%
  select(chrom, start, end, name = id, score, strand, thickStart, thickEnd, itemRgb) %>%
  arrange(chrom, start)
trackfile <- paste0(dir, "/analysis/assembly_dm6bal", assembly, ".bed")
cat(paste0('track name="Balancer assembly (dm6bal', assembly, ')" visibility=1 itemRgb=on\n'), file = trackfile)
write.table(as.track, file = trackfile, sep = "\t", quote = F, row.names = F, col.names = F, append = T)


calculate_overlap <- function(breakpoint, ref.coords)
{
  if (length(ref.coords) != 2)
    return(0L)
  else
    return(breakpoint[match("lower", ref.coords)] - breakpoint[match("higher", ref.coords)])
}

br <- rbind(
  as %>% mutate(breakpoint = ifelse(strand == "-", breakpoint2 - 1L, breakpoint2), ref.coords = ifelse(strand == "-", "higher", "lower")) %>%
    select(id = id2, chrom, breakpoint, ref.coords),
  as %>% mutate(breakpoint = ifelse(strand == "+", breakpoint1 - 1L, breakpoint1), ref.coords = ifelse(strand == "+", "higher", "lower")) %>%
    select(id = id1, chrom, breakpoint, ref.coords)
) %>%
  arrange(chrom, breakpoint, ref.coords) %>%
  group_by(id, chrom) %>%
  mutate(duplication_size = calculate_overlap(breakpoint, ref.coords)) %>%
  ungroup

write.table(br, file = paste0(dir, "/analysis/breakpoints_dm6_v", assembly, ".tab"), sep = "\t", quote = F, row.names = F)


bed <- rbind(
  as %>% mutate(
    breakpoint = ifelse(strand == "+", breakpoint1 - 1L, breakpoint1),
    start = ifelse(strand == "+", breakpoint, breakpoint - margin),
    id = paste0(id1, ifelse(strand == "+", "_B", "_A"))
  ) %>%
    select(chrom, start, id, strand),
  as %>% mutate(
    breakpoint = ifelse(strand == "-", breakpoint2 - 1L, breakpoint2),
    start = ifelse(strand == "-", breakpoint, breakpoint - margin),
    id = paste0(id2, ifelse(strand == "-", "_B", "_A"))
  ) %>%
    select(chrom, start, id, strand)
) %>%
  mutate(end = start + margin, score = margin, strand = "+") %>% # bal.strand (i.e. "+", not strand) taken
  select(chrom, start, end, id, score, strand) %>%
  arrange(chrom, start, end)

write.table(bed, file = paste0(dir, "/analysis/breakpoints_dm6_v", assembly, ".bed"), sep = "\t", quote = F, col.names = F, row.names = F)


asb <- as %>%
  group_by(chrom) %>%
  mutate(chrom_size = max(breakpoint1, breakpoint2)) %>%
  ungroup %>%
  mutate(length = abs(breakpoint2 - breakpoint1) + 1L) %>%
  mutate(bal.chrom = unique(chrom)[cumsum(grepl("start$|end$", id1))]) %>%
  group_by(bal.chrom) %>%
  mutate(bal.breakpoint1 = cumsum(length) - length + 1L, bal.breakpoint2 = cumsum(length), bal.strand = "+") %>%
  mutate(bal.chrom_size = max(bal.breakpoint1, bal.breakpoint2)) %>%
  ungroup


bal.br <- rbind(
  asb %>% mutate(chrom = bal.chrom, breakpoint = bal.breakpoint2, ref.coords = ifelse(strand == "-", "higher", "lower")) %>%
    select(id = id2, chrom, breakpoint, ref.coords),
  asb %>% mutate(chrom = bal.chrom, breakpoint = bal.breakpoint1 - 1L, ref.coords = ifelse(strand == "+", "higher", "lower")) %>%
    select(id = id1, chrom, breakpoint, ref.coords)
) %>%
  unique %>%
  separate(id, c("id_1", "id_2"), "_") %>%
  group_by(id_1, chrom, breakpoint) %>%
  summarize(id = paste0(first(id_1), "_", paste(id_2, ref.coords, sep = "_ref.", collapse = "_"))) %>%
  ungroup %>%
  select(id, chrom, breakpoint) %>%
  arrange(chrom, breakpoint)

write.table(bal.br, file = paste0(dir, "/analysis/breakpoints_dm6bal", assembly, ".tab"), sep = "\t", quote = F, row.names = F)


bal.bed <- rbind(
  asb %>% mutate(
    breakpoint = bal.breakpoint2,
    start = ifelse(bal.strand == "+", breakpoint - margin, breakpoint),
    id = paste0(id2, ifelse(strand == "-", "_B", "_A"))
  ) %>%
    select(bal.chrom, start, id, strand),
  asb %>% mutate(
    breakpoint = bal.breakpoint1 - 1L,
    start = ifelse(bal.strand == "+", breakpoint, breakpoint - margin),
    id = paste0(id1, ifelse(strand == "+", "_B", "_A"))
  ) %>%
    select(bal.chrom, start, id, strand)
) %>%
  mutate(end = start + margin, score = margin) %>% # strand (not bal.strand) taken
  select(chrom = bal.chrom, start, end, id, score, strand) %>%
  arrange(chrom, start, end)

write.table(bal.bed, file = paste0(dir, "/analysis/breakpoints_dm6bal", assembly, ".bed"), sep = "\t", quote = F, col.names = F, row.names = F)


dm6_chrom_sizes <- tbl_df(data.table::fread(paste0(dir, "/liftOver/dm6.chrom.sizes"),
  header = F, col.names = c("chrom", "chrom_size"), data.table = F))

dm6_unchanged_chain <- dm6_chrom_sizes %>%
  filter(!grepl("^chr[23]", chrom)) %>%
  mutate(chain = paste0(
    "chain 1000 ", # score -- chain score
    chrom, " ", # tName -- chromosome (reference sequence)
    chrom_size, " ", # tSize -- chromosome size (reference sequence)
    "+ 0 ", # tStrand -- strand (reference sequence), tStart -- alignment start position (reference sequence)
    chrom_size, " ", # tEnd -- alignment end position (reference sequence)
    chrom, " ", # qName -- chromosome (query sequence)
    chrom_size, " ", # qSize -- chromosome size (query sequence)
    "+ 0 ", # qStrand -- strand (query sequence), qStart -- alignment start position (query sequence)
    chrom_size, " ", # qEnd -- alignment end position (query sequence)
    100000 + row_number(), "\n", # id -- chain ID
    chrom_size, "\n" # size -- the size of the ungapped alignment
  ))

asb_chain <- asb %>%
  mutate(chain_from_dm6 = paste0(
    "chain 1000 ", # score -- chain score
    chrom, " ", # tName -- chromosome (reference sequence)
    chrom_size, " ", # tSize -- chromosome size (reference sequence)
    "+ ", # tStrand -- strand (reference sequence)
    ifelse(strand == "+", breakpoint1 - 1L, breakpoint2 - 1), " ", # tStart -- alignment start position (reference sequence)
    ifelse(strand == "+", breakpoint2, breakpoint1), " ", # tEnd -- alignment end position (reference sequence)
    bal.chrom, " ", # qName -- chromosome (query sequence)
    bal.chrom_size, " ", # qSize -- chromosome size (query sequence)
    strand, " ", # qStrand -- strand (query sequence)
    ifelse(strand == "+", bal.breakpoint1 - 1L, bal.chrom_size - bal.breakpoint2), " ", # qStart -- alignment start position (query sequence)
    ifelse(strand == "+", bal.breakpoint2, bal.chrom_size - bal.breakpoint1 + 1L), " ", # qEnd -- alignment end position (query sequence)
    row_number(), "\n", # id -- chain ID
    length, "\n" # size -- the size of the ungapped alignment
  ),
  chain_to_dm6 = paste0(
    "chain 1000 ", # score -- chain score
    bal.chrom, " ", # tName -- chromosome (reference sequence)
    bal.chrom_size, " ", # tSize -- chromosome size (reference sequence)
    "+ ", # tStrand -- strand (reference sequence)
    bal.breakpoint1 - 1L, " ", # tStart -- alignment start position (reference sequence)
    bal.breakpoint2, " ", # tEnd -- alignment end position (reference sequence)
    chrom, " ", # qName -- chromosome (query sequence)
    chrom_size, " ", # qSize -- chromosome size (query sequence)
    strand, " ", # qStrand -- strand (query sequence)
    ifelse(strand == "+", breakpoint1 - 1L, chrom_size - breakpoint1), " ", # qStart -- alignment start position (query sequence)
    ifelse(strand == "+", breakpoint2, chrom_size - breakpoint2 + 1L), " ", # qEnd -- alignment end position (query sequence)
    row_number(), "\n", # id -- chain ID
    length, "\n" # size -- the size of the ungapped alignment
  ))

writeLines(c(asb_chain$chain_from_dm6, dm6_unchanged_chain$chain), paste0(dir, "/liftOver/dm6ToDm6bal", assembly, ".over.chain"))
writeLines(c(asb_chain$chain_to_dm6, dm6_unchanged_chain$chain), paste0(dir, "/liftOver/dm6bal", assembly, "ToDm6.over.chain"))
system(paste0("gzip -f ", dir, "/liftOver/dm6ToDm6bal", assembly, ".over.chain ", dir, "/liftOver/dm6bal", assembly, "ToDm6.over.chain"))


chrom_sizes <- rbind(
  asb %>%
    select(chrom, chrom_size) %>%
    unique,
  dm6_unchanged_chain %>%
    select(chrom, chrom_size)
) %>%
  arrange(desc(chrom_size))

bal.chrom_sizes <- rbind(
  asb %>%
    select(chrom = bal.chrom, chrom_size = bal.chrom_size) %>%
    unique,
  dm6_unchanged_chain %>%
    select(chrom, chrom_size)
) %>%
  arrange(desc(chrom_size))

# check that dm6 chrom size was properly calculated
stopifnot(chrom_sizes == dm6_chrom_sizes)

write.table(bal.chrom_sizes, file = paste0(dir, "/liftOver/dm6bal", assembly, ".chrom.sizes"), sep = "\t", quote = F, col.names = F, row.names = F)
