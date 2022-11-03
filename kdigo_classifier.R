kdigo_classifier <- function(df, date_format="%Y-%m-%d", basal="min", week = FALSE){
    # kdigo_classifier uses creatinine data from hospitalized patients to infer 
    # acute kidney injury based on kdigo guidelines
    #
    # stage 0: no AKI
    # stage 1:1.5-2.0 baseline OR increase in 0.3mg/dl in 24/48h
    # stage 2: 2.0-3.0 baseline
    # stage 3: >3.0 baseline OR max value > 4 mg/dl
    # NA: no date available
    #
    # Input: data frame with alternate columns with exam data and creatinine value
    # Each row is a patient, each column, date/value
    # Example:  2020-02-14    1.16    2020-02-17    1.27
    #           2020-03-25    0.77    2020-03-26    1.92
    #           2020-03-29    0.85    2020-04-01    1.57
    # Output : A numeric vector of AKI stage, each element respective to a patient
    # Example: 3  3  3  2  0  3
    # Parameters:
    #       date_format:    The date format utilized in the df, with separators.
    #                       Standard: "%Y-%m-%d"
    #       basal:          Baseline is the "first" value or the "min" value?
    #                       Standard: "min"
    #       week:           Iterate in a time-spam of a week (TRUE) or consider
    #                           all the data (FALSE)
    # Usage:    classified <- kdigo_classifier(df, date_format="%Y-%m-%d", basal = "min", week = FALSE)
    
    require(tidyverse)
    #setting functions
    df_fix <- function(df, date_format = date_format){
        header_data <- paste("date", seq(1,ncol(df)/2), sep = "_")
        header_value <-  paste("value", seq(1,ncol(df)/2), sep = "_")
        header <- c(rbind(header_data,header_value))
        colnames(df) <- header
        
        #converting to date
        df[,seq(1, ncol(df), 2)] <- apply(df[,seq(1, ncol(df), 2)] , MARGIN = 2, FUN = as.Date, format = date_format)
        
        df$pacient = as.character(rownames(df))
        
        df2 <- pivot_longer(df, 
                            cols = !pacient, 
                            names_to = '.value', 
                            names_pattern = '(date|value)')
        df3 <- df2 %>%
            mutate(across(c(pacient, date, value), as.numeric))
        
        df3 <- df3 %>% 
            arrange(pacient, date)
        
        df3 <- df3 %>% 
            filter(!is.na(date)) %>% 
            filter(!is.na(value))
        
        return(df3)
    }
    
    get_diff_date <- function(df_fixed_date){
        x <- df_fixed_date
        x1 <- append(x[1], x)
        x2 <- append(x, x[length(x)])
        xdiff <- x2-x1
        xdiff <- xdiff[-c(length(xdiff))]
        return(xdiff)
    }
    
    get_status <- function(df_kdigo, basal = basal){
        test <- df_kdigo
        max_val <- max(test$value)
        min_val <- min(test$value)
        start_val <- test$value[1]
        end_val <- test$value[nrow(test)]
        
        #what is basal
        if(basal == "min"){
            basal <- min_val
        } else if(basal == "first"){
            basal <- start_val
        } else{
            stop("Basal must be \"first\" or \"min\"")
        }
        
        
        if(max_val >= 4){ #set status as 3 by max val
            status <- status_final <- 3
        } else {        #set status as no lesion
            status <- status_final <- 0
        }
        
        for(i in seq(1,nrow(test))){
            if(i > 1){
                #status 1 if value >= 0.3 in 24/48h
                if(test$diff[i] <= 2){
                    if((test$value[i] - test$value[i-1]) >= 0.3){
                        status <- 1
                    } else if(i > 2){
                        if((test$value[i] - test$value[i-2]) >= 0.3){
                            status <- 1
                        }
                    }
                }
                if(isTRUE(week)){
                    if(val >= 4){ #status 3 val >= 3*basal
                        status <- 3
                    } else if(val >= 1.5*basal & val < 2*basal){ #status 1 val >= 1.5*basal & val < 2*basal in a week
                        status <- 1
                    } else if(val >= 2*basal & val < 3*basal){ #status 2 val >= 2*basal & val < 3*basal in a week
                        status <- 2
                    } else if(val >= 3*basal){ #status 2 val >= 2*basal & val < 3*basal in a week
                        status <- 3
                    }
                } else if(isFALSE(week)){
                    # get week interval
                    week <- FALSE
                    j <- i-1
                    while(isFALSE(week)){
                        if(positive_i(j) == 0){
                            week <- TRUE
                        }else if(sum(test$diff[seq(j,i)]) < 7){
                            j <- j - 1
                        }else if(sum(test$diff[seq(j,i)]) >= 7){
                            week <- TRUE
                            val<- max(test$value[seq(j,i)])
                            #status 3 val >= 4 & status == "1"
                            if(val >= 4){ #status 3 val >= 3*basal
                                status <- 3
                            } else if(val >= 1.5*basal & val < 2*basal){ #status 1 val >= 1.5*basal & val < 2*basal in a week
                                status <- 1
                            } else if(val >= 2*basal & val < 3*basal){ #status 2 val >= 2*basal & val < 3*basal in a week
                                status <- 2
                            } else if(val >= 3*basal){ #status 2 val >= 2*basal & val < 3*basal in a week
                                status <- 3
                            }
                        }
                    }
                }else{
                    error("Week should be TRUE or FALSE. Standard is FALSE.")
                }
                #based on basal
                val<- max(test$value[i])

            }
            if (status > status_final){
                status_final <- status
            }
        }
        
        return(status_final)
    }
    
    positive_i <- function(i){
        if(i > 0){
            return(i)
        } else {
            return(0)
        }
    }
    
    #setting code
    
    
    df_fixed <- df_fix(df)
    kdigo <- df_fixed %>%
        group_by(pacient) %>% 
        mutate(diff = get_diff_date(date)) %>% 
        mutate(stat = get_status(cur_data(), basal = basal)) %>% 
        distinct(pacient, .keep_all = TRUE)
    
    df_base <- data.frame(pacient = as.numeric(rownames(df)))
    final <- full_join(df_base, kdigo, by = "pacient")
    final_status <- final$stat
    return(final_status)
}
