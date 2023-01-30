import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


def create_label_lists(all_cross_test_label, all_cross_predictions):
    """
    :params all_cross_predictions: list of lists, in each list there are predictions (probabillities)
    :params all_cross_test_label: list of lists , in each list the labels of the according queries (hot one encoded as [Naa,2]
    :return tuple (label0_list,label1_list)
    """
    label0_list = []
    label1_list = []
    label0 = np.array([1, 0])
    label1 = np.array([0, 1])
    # print(label0)
    # print(all_cross_test_label[0][0])
    # print(np.all(label0,all_cross_test_label[0][0]))
    for i in range(len(all_cross_test_label)):
        for j in range(len(all_cross_test_label[i])):
            label = all_cross_test_label[i][j]
            if np.all(label == label0):  # label = 0
                label0_list.append(all_cross_predictions[i][j])
            elif np.all(label == label1):  # label = 1
                label1_list.append(all_cross_predictions[i][j])
            else:  # non valid query
                continue
    return label0_list, label1_list


def create_labels_histogram(label0_list, label1_list):
    """
    :params label0_list: list of the predictions with label = 0
    :params label1_list: list of the predictions with label = 1
    The function plots histogram of the labels distribution
    """
    label0_list = np.array(label0_list)
    label1_list = np.array(label1_list)
    sns.kdeplot(data=label0_list, label="non UBD's")
    sns.kdeplot(data=label1_list, label="UBD's")
    plt.xlabel("binding probability")
    plt.legend()
    plt.show()


# def find_percentile(percentile, all_cross_predictions):
#     """
#     :params perenctile
#     :params data - numpy list of numpy lists
#     """
#     new_arr = []
#     for i in range(len(all_cross_predictions)):
#         for j in range(len(all_cross_predictions[i])):
#             new_arr.append(all_cross_predictions[i][j])
#     predictions = np.concatenate(all_cross_predictions)
#     return np.percentile(predictions, percentile)
#

def find_threshold(precision, label0_list, label1_list, all_predictions):
    """
    :params: precision: value of precision to make
    :params label0_list: list of the predictions with label = 0
    :params label1_list: list of the predictions with label = 1
    The function returning the threshold needed to achieve precision
    """
    # all_predictions = np.concatenate(np.array([np.array(label0_list), np.array(label1_list)]))
    # all_predictions = np.sort(all_predictions)
    print(len(all_predictions))
    for i in range(len(all_predictions)):
        threshold = all_predictions[i]
        mask0 = label0_list > threshold
        num_0_elements = np.sum(mask0)
        mask1 = label1_list > threshold
        num_1_elements = np.sum(mask1)
        calculated_precision = (num_1_elements) / (num_0_elements + num_1_elements)
        if calculated_precision > precision:
            print("num1_elements :", num_1_elements)
            print("num0_elements :", num_0_elements)
            return threshold


# list1 = [[0, 1], [0, 1], [1, 0], [0, 1], [0, 0], [1, 0]]
# l1 = np.array([10, 7, 4])
# l2 = np.array([3, 2, 1])
# l3 = np.array([l1, l2])
# print(l3)
# #
# #
# print(find_percentile(33, l3))
# list1 = np.asarray(list1)
# list1 = np.array([list1])
# list2 = [0.5, 0.4, 0.3, 0.4, 0.9, 0.1]
# list2 = np.asarray(list2)
# list2 = np.array([list2])
#
# print(list1, list2)
#
# label0_list, label1_list = create_label_lists(list1, list2)
# create_labels_histogram(label0_list, label1_list)
