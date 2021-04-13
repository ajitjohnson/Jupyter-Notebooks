#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 20 20:39:36 2021
@author: Ajit Johnson Nirmal
Twitter Post Finder
"""

# Lib
import tweepy
import datetime
import pandas as pd


# Authentication
auth = tweepy.AppAuthHandler('VhVYh9RVGivvFanwDr7kg8qBt', 'oXG5rJDopPVgThiq6xEU3fKNonpkDTFMvazQjSvgkw66NeuVsf')

# test
api = tweepy.API(auth)
    
# Run Function
pages = ['NewsfromScience','ScienceNews','newscientist','Immunoscope',
         'nature','NatureMedicine','nresearchnews','natBME','NatureBiotech','NatureComms','NatImmunol','NatureNews','naturemethods',
         'NatureHumBehav','NatureAstronomy','NatureGenet','NatureNV','NatureClimate','NatureCellBio',
         'CellCellPress','TrendsCellBio','Dev_Cell','CellSystemsCP','CellPressNews','TrendsMolecMed','TrendsImmuno',
         'MolecularCell','CellStemCell','NeuroCellPress','Cancer_Cell','cellhostmicrobe','JCellBiol',
         'NewsfromScience','SciImmunology','ScienceMagazine','ScienceAdvances','scisignal','ScienceTM']
startDate = datetime.datetime(2021, 3, 21, 0, 0, 0) # Y | M | D | hour | Min | Second
endDate =   datetime.datetime.now()

output = go_viral(profile_ids=pages, post_count=100, startDate=startDate, endDate=endDate)
suffix = str(endDate.day) + '_' + str(endDate.month) + '_' + str(endDate.year) + '.csv'
output.to_csv('/Users/aj/Dropbox (Partners HealthCare)/Papers to read/twitter/' + suffix, index=False)

    
def go_viral (profile_ids, post_count=100, startDate=None, endDate=None):
    
    # Make a mega list of the posts
    all_posts = [] # intialize
    for pages in profile_ids:
        try:
            single_page = api.user_timeline(pages, count=post_count, tweet_mode="extended")
            print ('Processing: ' + pages)
            all_posts.extend(single_page)
        except Exception:
            print ('Failed Processing: ' + pages)
            pass
    
    # remove tweets that does not fall into the user defined time range
    if startDate is not None or endDate is not None:
        print ('Processing tweets within date range')
        
    if startDate is not None:
        tweets = []
        for tweet in all_posts:
            if tweet.created_at > startDate:
                tweets.append(tweet)

    if endDate is not None:
        tweets = []
        for tweet in all_posts:
            if tweet.created_at < endDate:
                tweets.append(tweet)

    if startDate is not None and endDate is not None:
        tweets = []
        for tweet in all_posts:
            if tweet.created_at > startDate and tweet.created_at < endDate:
                tweets.append(tweet)

    if startDate is None and endDate is None:
        tweets = all_posts
    
    data = [] # intialize to hold data
    
    print('Consolidating results')
    
    # create a dataframe with likes retweet page tweetid and date
    feature_list = ['id', 'created_at','retweet_count', 'favorite_count', 'full_text']   
    for features in feature_list:
         foi = 'status.'+ str(features)        
         tmp_feature = []
         for status in tweets:
             tmp_feature.append(eval(foi))
         # add the tmp feature to the main feature list
         data.append(tmp_feature)
    
    # extracting user data
    user_feature_list = [ "name",  "screen_name",  "followers_count"]
    for features in user_feature_list:
         foi = 'status.user.'+ str(features)        
         tmp_feature = []
         for status in tweets:
             tmp_feature.append(eval(foi))
         # add the tmp feature to the main feature list
         data.append(tmp_feature)

    # convert to dataframe
    combined_data = pd.DataFrame(data, index=feature_list + user_feature_list).T
    
    # Make sure there is atleast one like and retweet
    combined_data[['retweet_count', 'favorite_count']] = combined_data[['retweet_count', 'favorite_count']].replace(0, 1)
    
    # remove re-tweets
    to_remove = combined_data[combined_data['full_text'].str.startswith('RT')].index
    combined_data = combined_data.drop(to_remove)
        
    # Engagement score
    combined_data['raw_ES'] = (combined_data['retweet_count'] * 2) * combined_data['favorite_count']
    combined_data['norm_ES'] = (combined_data['raw_ES'] / combined_data['followers_count']) * 100
    combined_data['Percent Engagement'] = ((combined_data['retweet_count'] + combined_data['favorite_count']) / combined_data['followers_count']) * 100
    
    # AJ Enrichment Score
    # Account for profile specific engagement 
    retweet_sum = combined_data.groupby('screen_name')['retweet_count'].transform('sum')
    favorite_sum = combined_data.groupby('screen_name')['favorite_count'].transform('sum')
    
    # normalize the re tweets and likes
    re_tweets_n = combined_data['retweet_count'].div(retweet_sum)
    favorite_n = combined_data['favorite_count'].div(favorite_sum)
    
    # generate a score 
    score = ((re_tweets_n * 2) + favorite_n) * 100
    
    # add the score t0 the dataframe
    combined_data ['aj_score'] = score
    
    # add a link column
    combined_data ['link'] = 'https://twitter.com/twitter/status/' + combined_data['id'].astype(str)
    
    # reorder the columns
    re_column = ["screen_name",'link','retweet_count', 'favorite_count', 'full_text', 'norm_ES','Percent Engagement','aj_score']
    combined_data = combined_data[re_column]
    
    # sort
    combined_data = combined_data.sort_values(by='aj_score', ascending=False)
    
    # return 
    return combined_data
    
    
    
    
    
    
    
    
    
    
