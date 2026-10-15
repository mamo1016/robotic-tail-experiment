function [eeg_feature, foot_feature] = feature_ext_sig(all_sub_diffs)

sub_list = fieldnames(all_sub_diffs);
sub_size = length(sub_list);
eeg_feature = zeros(sub_size,1200);
foot_feature = zeros(sub_size,1200);
for sub = 1:sub_size    
    eeg_feature(sub,1:400)      = all_sub_diffs.(sub_list{sub}).eeg_significant.expected_vs_expected_duplicated(:);
    eeg_feature(sub,401:800)    = all_sub_diffs.(sub_list{sub}).eeg_significant.unexpected_vs_expected(:);
    eeg_feature(sub,801:1200)   = all_sub_diffs.(sub_list{sub}).eeg_significant.right_after_vs_right_before(:);

    foot_feature(sub,1:400)      = all_sub_diffs.(sub_list{sub}).foot_significant.expected_vs_expected_duplicated(:);
    foot_feature(sub,401:800)    = all_sub_diffs.(sub_list{sub}).foot_significant.unexpected_vs_expected(:);
    foot_feature(sub,801:1200)   = all_sub_diffs.(sub_list{sub}).foot_significant.right_after_vs_right_before(:);
end