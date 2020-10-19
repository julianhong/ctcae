#process.R
#translates cTakes SNOMED output into CTCAE
#based on manual review of 100 rad onc notes
#copyright Duke University, Julian Hong
#GNU GPL v2.0 license


#load notes extracts

load("<file path>")


#import the vocabulary file from OHDSI Athena
concept_relationships <- read.table("<file path>", sep = "\t", dec = ".", header = TRUE,
                                    encoding = "UTF-8", stringsAsFactors = FALSE, quote = "")

#concept file needed to have "#" removed for reading
concepts <- read.table("<file path>", sep = "\t", dec = ".", header = TRUE,
                                    encoding = "UTF-8", stringsAsFactors = FALSE, quote = "")

#now separate out the SNOMED and MedDRA codes
library(dplyr) #v0.8.5

snomed <- filter(concepts, vocabulary_id == "SNOMED")
meddra <- filter(concepts, vocabulary_id == "MedDRA")

#Generate the map
map <- filter(concept_relationships, relationship_id == "SNOMED - MedDRA eq")
temp <- map %>%
  left_join( snomed, by = c("concept_id_1" = "concept_id")) %>%
  left_join(meddra, by = c("concept_id_2" = "concept_id"))

#verified join and cleanup
map <- temp
rm(temp)

#convert to characters
map$concept_code.x <- as.character(map$concept_code.x)
map$concept_code.y <- as.character(map$concept_code.y)

#now let's join snomed to meddra
temp <- left_join(notes, map, by = c("code" = "concept_code.x"))

#convert this to a long form by individual concepts
library(tidyr) #v1.0.0
ctakeslong <-
  temp %>%
  select(NOTE_ID, concept_code.y, polarity, preferredText, textsem, concept_name.y) %>%
  group_by(NOTE_ID, concept_code.y) %>%
  summarize(polarity = max(as.numeric(as.character(polarity))), firstsnomed = first(preferredText), firstsem = first(textsem), meddra = first(concept_name.y))

#there were no discordances within single notes

rm(temp)

#let's now take a look at concepts that were not mapped to Meddra to identify potential missed symptoms (none were symptoms)
ctakesfilt <-
  ungroup(ctakeslong) %>%
  filter(is.na(concept_code.y)) %>%
  select(firstsnomed, firstsem) %>%
  group_by(firstsnomed) %>%
  summarize(sem = first(firstsem), num = n())

summary(as.factor(ctakesfilt$sem))

#import reviewer data
revA <- read.csv("<revA file path>", check.names = FALSE, stringsAsFactors = FALSE)
revA$NOTE_ID <- as.factor(revA$NOTE_ID)
revA$reviewer <- "revA"

revB <- read.csv("<revB file path>", check.names = FALSE, stringsAsFactors = FALSE)
revB$NOTE_ID <- as.factor(revB$NOTE_ID)
revB$reviewer <- "revB"

#check the NOTE_IDs
revA$NOTE_ID == revB$NOTE_ID

#let's clean up here and get rid of the reviwer column
revA <- select(revA, -reviewer)
revB <- select(revB, -reviewer)

#convert all the NAs to 0
revA[is.na(revA)] <- 0
revB[is.na(revB)] <- 0


#alright let's import our consensus data now
consensus <- read.csv("<file path>", check.names = FALSE, stringsAsFactors = FALSE)

#similar to the ctakes data we're going to have to fill out all the blanks here too to match the columns/rows

#first pull column names and verify our import
colnames(revB) == colnames(revA) #confirm our column names are the same
revB$NOTE_ID == revA$NOTE_ID #confirms our rows are the same
revBcol <- colnames(revB) #we'll just use the column names to conform the consensus dataframe which was simplified for review

#so what are we missing in the consensus frame (which were all non-present symptoms by both reviewers)
missing <- setdiff(revBcol, colnames(consensus)) #we want to find the revB columns that are missing for pertinent ruleout

#now we have to readd in the columns that were in the original dataframe though
#in this case all of the columns are a subset of the full thing
#make a new 0 matrix to add first
add <- matrix(0L, nrow = nrow(consensus), ncol = length(missing))
colnames(add) <- missing
temp1 <- cbind(as.data.frame(consensus), add, stringsAsFactors = FALSE)
consensus <- rbind(revB[0,], temp1[,revBcol]) #combine the column header here with the new stuff

#verify comparison
setdiff(revBcol, colnames(consensus))
sum(match(revBcol, colnames(consensus)) == 1:ncol(consensus)) #should be 838

#we have to do the same thing for rows which should be a left join to resort
consensus$NOTE_ID <- as.factor(consensus$NOTE_ID)
consensus <- left_join(select(revB, NOTE_ID), consensus)

#fill those in with 0s now
consensus[is.na(consensus)] <- 0

rm(missing)

###
#There are relevant extracted MEDDRA codes that weren't in the CTCAE term dictionary. Let's take a look at those and get a plaintext explanation for them
#there are 145 codes here
#CTCAE is from the NCI https://ctep.cancer.gov/protocolDevelopment/electronic_applications/ctc.htm (v5)

ctakescol <- unique(na.omit(ctakeslong$concept_code.y))

ctakesnotctcae <-
  as.data.frame(setdiff(ctakescol, revBcol)) %>%
  left_join(map, by = c("setdiff(ctakescol, revBcol)" = "concept_code.y"))

#this will require a ctakes dictionary which we can build off of ctakeslong (which contains all of the ctakes codes mapped to Meddra)
#let's compile a list to see if their synonyms are mapped already
ctakesterms <-
  ctakeslong %>%
  filter(!is.na(concept_code.y)) %>%
  group_by(concept_code.y) %>%
  summarize(firstsem = first(firstsem), meddra = first(meddra), num = n())

ctakesterms$ctcae[ctakesterms$concept_code.y %in% revBcol] <- 1 #1 if yes, 0 if not in CTCAE
ctakesterms$ctcae[!(ctakesterms$concept_code.y %in% revBcol)] <- 0 #1 if yes, 0 if not in CTCAE

sum(ctakesterms$ctcae) #output out the # in CTCAE for reference

#write it out
write.csv(ctakesterms, paste("<destination directory>", date, "/ctakesterms.csv", sep = ""))

#JH has prepared a thesaurus to convert these terms and consolidate with CTCAE
thesaurus <- read.csv("<thesaurus file>", check.names = FALSE, stringsAsFactors = FALSE, colClasses = rep("character", 5))

#let's make some changes to ctakes long
#we want to basically convert the existing codes if they are in our thesaurus
ctakeslongsyn <- left_join(ctakeslong, thesaurus, by = c("concept_code.y" = "concept_code"))

#anything with a synonym should get written into concept_code.y
#create a ctakes long analog which has the new features; should be the same length as all detected concept-codes
ctakeslongsyn$concept_code.y[!is.na(ctakeslongsyn$synonymcode)] <- ctakeslongsyn$synonymcode[!is.na(ctakeslongsyn$synonymcode)]

ctakeslongsyn <- #this part will shorten us down to < original by collapsing codes
  ctakeslongsyn %>%
  select(colnames(ctakeslong)) %>%
  group_by(NOTE_ID, concept_code.y) %>%
  summarize(polarity = max(as.numeric(as.character(polarity))), firstsnomed = first(firstsnomed), firstsem = first(firstsem), meddra = first(meddra))

ctakes <-
  select(ctakeslongsyn, NOTE_ID, concept_code.y, polarity) %>%
  spread(concept_code.y, polarity)

ctakescol <- colnames(ctakes)#reset the column names here goes from 182-> 174

#now we need to add in all the other missing symptoms that were not detected though via cTakes
#it's also OK to eliminate those codes without CTCAE mapping at this point but let's take a look and see what they are now.
#Make sure any meaningful symptoms with > 1 mention are coded

temp <-
  ctakeslongsyn %>%
  filter(!is.na(concept_code.y)) %>%
  group_by(concept_code.y) %>%
  summarize(firstsem = first(firstsem), meddra = first(meddra), num = n())
temp$ctcae[temp$concept_code.y %in% revBcol] <- 1 #1 if yes, 0 if not in CTCAE
temp$ctcae[!(temp$concept_code.y %in% revBcol)] <- 0 #1 if yes, 0 if not in CTCAE

write.csv(temp, paste("<output file directory path>", date, "/ctakestermspostmap.csv", sep = ""))
rm(temp)

#so what are we missing
missing <- setdiff(revBcol, ctakescol) #we want to find the revB columns that are missing for pertinent ruleout

#now we want to set revB's headers
temp <-ctakes[, intersect(revBcol, ctakescol)]

#now we have to readd in the columns that were in the original dataframe though
#make a new 0 matrix
add <- matrix(0L, nrow = nrow(temp), ncol = length(missing))
colnames(add) <- missing
temp1 <- cbind(as.data.frame(temp), add, stringsAsFactors = FALSE)
ctakes <- rbind(revB[0,], temp1[,revBcol]) #combine the column header here with the new stuff

#verify comparison
setdiff(revBcol, colnames(ctakes))
sum(match(revBcol, colnames(ctakes)) == 1:ncol(ctakes)) #should be 838

#we have to do the same thing for rows which should be easy with a left join to resort
ctakes$NOTE_ID <- as.factor(ctakes$NOTE_ID)
ctakes <- left_join(select(revB, NOTE_ID), ctakes)

#fill those in with 0s now
ctakes[is.na(ctakes)] <- 0

#write out the ctakes output file
write.csv(ctakes, paste("<out file directory path>", date, "/ctakes.csv", sep = ""))

###
#now that all our data (revA, revB, consensus, and cTakes) is initialized, we also need to do some grouping based on expert discrepancies
#(Fairchild, IJROBP 2020)
#we have a grouping dictionary, but we'll just do this functionally here.

#initialize known synonyms (8 categories with 37; net -29 variables) -> 809 variables
#let's function this

groupsx <- function(ssxtab) {
  ssxtabtrans <- ssxtab #initialize

  #hematochezia (2)
  #10009887 (colitis - revA)
  #10038064 (rectal hemorrhage - revB)
  ssxtabtrans$hematochezia[ssxtab$`10009887`== -1 | ssxtab$`10038064` == -1] <- -1
  ssxtabtrans$hematochezia[ssxtab$`10009887`==1 | ssxtab$`10038064` == 1] <- 1 #will finish with positives since if both mentioned that is default (just for safety)

  #now drop
  ssxtabtrans <- select(ssxtabtrans, -c(`10009887`, `10038064`))

  #vision changes (2)
  #10005886 blurred vision
  #10047516 visual acuity

  ssxtabtrans$vision[ssxtab$`10005886` == -1 | ssxtab$`10047516` == -1] <- -1
  ssxtabtrans$vision[ssxtab$`10005886` == 1 | ssxtab$`10047516` == 1] <- 1

  ssxtabtrans <- select(ssxtabtrans, -c(`10005886`, `10047516`))


  #pain 10033371 (20)
  #10000081 abdominal pain
  #10002167 anal pain
  #10015388 esophageal pain
  #10024561 lip pain

  #10054520 oral dysesthesia
  #10031009 oral pain
  #10044055 tooth pain

  #10016059 face pain
  #10062501 non-cardiac chest pain
  #10008496 CW pain
  #10035623 pleuritic chest pain
  #10003988 back pain
  #10028836 neck pain
  #10033425 extremity pain
  #10006298 breast pain

  #10068319 throat pain
  #10041367 sore throat

  #10033474 skin pain
  #10049120 scalp pain
  ssxtabtrans$pain[ssxtab$`10033371` == -1 | ssxtab$`10000081` == -1 | ssxtab$`10002167` == -1 | ssxtab$`10015388` == -1 | ssxtab$`10024561` == -1 | ssxtab$`10054520` == -1 |
                     ssxtab$`10031009` == -1 | ssxtab$`10044055` == -1 | ssxtab$`10016059` == -1 | ssxtab$`10062501` == -1 | ssxtab$`10008496` == -1 |
                     ssxtab$`10035623` == -1 | ssxtab$`10003988` == -1 | ssxtab$`10028836` == -1 | ssxtab$`10033425` == -1 | ssxtab$`10006298` == -1 |
                     ssxtab$`10068319` == -1 | ssxtab$`10041367` == -1 | ssxtab$`10033474` == -1 | ssxtab$`10049120` == -1] <- -1

  ssxtabtrans$pain[ssxtab$`10033371` == 1 | ssxtab$`10000081` == 1 | ssxtab$`10002167` == 1 | ssxtab$`10015388` == 1 | ssxtab$`10024561` == 1 | ssxtab$`10054520` == 1 |
                     ssxtab$`10031009` == 1 | ssxtab$`10044055` == 1 | ssxtab$`10016059` == 1 | ssxtab$`10062501` == 1 | ssxtab$`10008496` == 1 |
                     ssxtab$`10035623` == 1 | ssxtab$`10003988` == 1 | ssxtab$`10028836` == 1 | ssxtab$`10033425` == 1 | ssxtab$`10006298` == 1 |
                     ssxtab$`10068319` == 1 | ssxtab$`10041367` == 1 | ssxtab$`10033474` == 1 | ssxtab$`10049120` == 1] <- 1

  ssxtabtrans <- select(ssxtabtrans, -c(`10033371`, `10000081`, `10002167`, `10015388`, `10024561`, `10054520`, `10031009`, `10044055`, `10016059`, `10062501`, `10008496`,
                                        `10035623`, `10003988`, `10028836`, `10033425`, `10006298`, `10068319`, `10041367`, `10033474`, `10049120`))

  #weakness (5)
  #10062572 generalized muscle weakness
  #10065776 weakness lower limb
  #10065895 weakness upper limb
  #10065780 weakness left
  #10065794 weakness right
  ssxtabtrans$weak[ssxtab$`10062572` == -1 | ssxtab$`10065776` == -1 | ssxtab$`10065895` == -1 | ssxtab$`10065780` == -1 | ssxtab$`10065794` == -1] <- -1
  ssxtabtrans$weak[ssxtab$`10062572` == 1 | ssxtab$`10065776` == 1 | ssxtab$`10065895` == 1 | ssxtab$`10065780` == 1 | ssxtab$`10065794` == 1] <- 1

  ssxtabtrans <- select(ssxtabtrans, -c(`10062572`, `10065776`, `10065895`, `10065780`, `10065794`))

  #edema (2)
  #10050068 limb edema
  #10025233 lymphedema
  ssxtabtrans$edema[ssxtab$`10050068` == -1 | ssxtab$`10025233` == -1] <- -1
  ssxtabtrans$edema[ssxtab$`10050068` == 1 | ssxtab$`10025233` == 1] <- 1

  ssxtabtrans <- select(ssxtabtrans, -c(`10050068`, `10025233`))

  #altered mental status (2)
  #10009845 cognitive disturbance
  #10012373 depressed consciousness
  ssxtabtrans$ams[ssxtab$`10009845` == -1 | ssxtab$`10012373` == -1] <- -1
  ssxtabtrans$ams[ssxtab$`10009845` == 1 | ssxtab$`10012373` == 1] <- 1

  ssxtabtrans <- select(ssxtabtrans, -c(`10009845`, `10012373`))

  #facial weakness (2)
  #10051272 facial muscle weakness
  #10061457 facial nerve disorder
  ssxtabtrans$faceweak[ssxtab$`10051272` == -1 | ssxtab$`10061457` == -1] <- -1
  ssxtabtrans$faceweak[ssxtab$`10051272` == 1 | ssxtab$`10061457` == 1] <- 1

  ssxtabtrans <- select(ssxtabtrans, -c(`10051272`, `10061457`))

  #cough (2)
  #10011224 cough
  #10036790 productive cough
  ssxtabtrans$cough[ssxtab$`10011224` == -1 | ssxtab$`10036790` == -1] <- -1
  ssxtabtrans$cough[ssxtab$`10011224` == 1 | ssxtab$`10036790` == 1] <- 1

  ssxtabtrans <- select(ssxtabtrans, -c(`10011224`, `10036790`))

  #might need some kind of combination for skin reaction

  ssxtabtrans[is.na(ssxtabtrans)] <- 0

  return(ssxtabtrans)
}

#Final output here:
revAtrans <- groupsx(revA)
revBtrans <- groupsx(revB)
ctakestrans <- groupsx(ctakes)
consensustrans <- groupsx(consensus)
